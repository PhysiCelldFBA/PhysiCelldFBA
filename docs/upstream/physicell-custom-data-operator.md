# Proposed PhysiCell fix: `Custom_Cell_Data::operator[]` silently corrupts data

Draft for a pull request / issue against upstream PhysiCell. Written against **PhysiCell 1.14.2**.

---

## Summary

Writing to a `custom_data` name that is not declared in the cell definition silently overwrites
the **first declared custom variable** instead of reporting an error. There is no warning, no
exception and no non-zero exit code. The corrupted variable keeps being overwritten on every
write, so a single typo in a user's `custom.cpp` produces plausible-looking but wrong output for
the entire run.

A secondary effect makes it worse: the failed lookup *inserts* the bad name into the lookup map,
after which `find_variable_index()` returns `0` for that name rather than `-1` — defeating the one
function a defensive user would reach for.

---

## Location

`core/PhysiCell_custom.cpp`, lines 236-239:

```cpp
double& Custom_Cell_Data::operator[]( std::string name )
{
	return variables[ name_to_index_map[name] ].value; 
}
```

with, in `core/PhysiCell_custom.h` line 110:

```cpp
std::unordered_map<std::string,int> name_to_index_map; 
```

## Mechanism

`name_to_index_map` is a `std::unordered_map`. Calling `operator[]` on an associative container
with a key that is not present **default-constructs and inserts** the value — which for `int` is
`0`. So for any undeclared `name`:

1. `name_to_index_map[name]` inserts `{name, 0}` and returns `0`;
2. the function returns a reference to `variables[0].value`;
3. the caller's write lands on whichever variable happens to be declared first.

The adjacent function already anticipates this class of problem and even documents the intent in a
comment (lines 201-208):

```cpp
int Custom_Cell_Data::find_variable_index( std::string name )
{
	// this should return -1 if not found, not zero 
	auto out = name_to_index_map.find( name ); 
	if( out != name_to_index_map.end() )
	{ return out->second; }
	return -1; 
}
```

`find_variable_index` is correct in isolation. But because `operator[]` has already inserted the
bad key, a subsequent `find_variable_index(name)` **finds** it and returns `0`. The guard a user
writes to protect themselves therefore reports that the variable exists, at index 0.

`find_vector_variable_index` is unaffected: it uses a linear scan over `vector_variables` (the
map-based version is commented out at lines 210-215).

## Reproduction

```cpp
// cell definition declares, in this order: "my_important_value", "other"
pCell->custom_data["my_important_value"] = 1.0;
pCell->custom_data["typo_in_this_name"]  = 42.0;   // undeclared

std::cout << pCell->custom_data["my_important_value"];        // prints 42
std::cout << pCell->custom_data.find_variable_index("typo_in_this_name"); // prints 0, not -1
```

Expected: an error on the second write, or at minimum a warning and no side effect.
Actual: `my_important_value` is destroyed, and the subsequent index lookup lies.

**Verified on 1.14.2**, compiling `core/PhysiCell_custom.cpp` standalone against a minimal driver:

```
before:  first=1 second=2
find_variable_index("nope") = -1   (expect -1)
after writing to undeclared "nope":
  first  = 42   (expect 1, bug gives 42)
  second = 2    (expect 2)
  find_variable_index("nope") = 0  (expect -1, bug gives 0)
```

## Impact

- Silent, total corruption of the first declared custom variable, for the whole run.
- The failure is invisible in output: the variable contains well-formed numbers of the right
  magnitude, just belonging to a different quantity.
- The documented defensive check (`find_variable_index`) stops working after the first bad write.
- A long-lived workaround in affected projects is to declare a throwaway variable first so it
  absorbs stray writes. That treats the symptom and hides the typo.
- `operator[]` mutates the container, so a function that reads like an accessor is in fact a
  writer. If a shared `Custom_Cell_Data` (for example on a `Cell_Definition`) is ever read
  concurrently, that is also a latent thread-safety issue.

In our own codebase this manifested as a project writing an exchange flux under the name
`R_EX_lac__D_e` while the configuration declared `R_EX_lac__L_e`. Every cell wrote its lactate
flux into the throwaway variable at every timestep, and the declared variable stayed at `0.0` in
every output file, for a long time, unnoticed.

## Proposed fix

Fail loudly, matching how PhysiCell already handles an unknown substrate name in
`core/PhysiCell_cell.cpp` (lines 2776, 2788, 2919 and others print an explanatory message and
`exit(-1)`):

```cpp
double& Custom_Cell_Data::operator[]( std::string name )
{
	auto search = name_to_index_map.find( name ); 
	if( search == name_to_index_map.end() )
	{
		std::cout << std::endl 
			<< "Error: custom data variable \"" << name << "\" does not exist." << std::endl
			<< "       Declare it in the <custom_data> block of the cell definition," << std::endl
			<< "       or test with find_variable_index() before accessing it." << std::endl
			<< "       Please check for typos!" << std::endl << std::endl; 
		exit(-1); 
	}
	return variables[ search->second ].value; 
}
```

The change is self-contained: one function, no signature change, no new dependency.

**Verified** by applying this patch to a copy of `core/PhysiCell_custom.cpp` and rebuilding the
driver above: correct accesses are unaffected (`first=7 second=8`, exit 0), the undeclared access
aborts with the message and a non-zero exit code, and `find_variable_index` keeps returning `-1`
because the map is never polluted.

Because the function returns `double&` it cannot report failure through its return value, so the
options are to abort, to throw, or to return a scratch reference. Aborting is proposed because it
matches the existing convention for configuration typos elsewhere in PhysiCell, and because this
is a programming error rather than a runtime condition.

### Less disruptive alternative

If maintainers prefer not to turn an existing (if broken) call into a fatal error, the same
corruption can be avoided non-fatally by warning once and directing the write to a scratch value:

```cpp
double& Custom_Cell_Data::operator[]( std::string name )
{
	auto search = name_to_index_map.find( name ); 
	if( search == name_to_index_map.end() )
	{
		static double scratch = 0.0; 
		static std::set<std::string> reported; 
		#pragma omp critical
		{
			if( reported.insert(name).second )
			{ 
				std::cout << "Warning: custom data variable \"" << name 
					<< "\" does not exist; the access is being ignored." << std::endl; 
			}
		}
		scratch = 0.0; 
		return scratch; 
	}
	return variables[ search->second ].value; 
}
```

This loses the write rather than misplacing it, and keeps existing runs alive. Note the scratch
reference is shared between threads, so it is only safe as a sink for discarded writes, never as
a value anyone reads.

### Compatibility

Any code that currently relies on the present behaviour is, by definition, writing to an
undeclared name and corrupting `variables[0]`. We are not aware of a legitimate use. The fatal
version will stop runs that were previously producing wrong numbers silently, which is the
intended outcome, but it is a behaviour change worth calling out in release notes.

Neither variant changes the signature, the header, or the behaviour of any correct call.

## Suggested regression test

```cpp
Custom_Cell_Data cd; 
cd.add_variable( "first" , "none" , 1.0 ); 
cd.add_variable( "second" , "none" , 2.0 ); 

// an undeclared name must not resolve to index 0 
assert( cd.find_variable_index( "nope" ) == -1 ); 

// ... and must still report -1 after an attempted access,
// i.e. the lookup map must not have been polluted 
assert( cd.find_variable_index( "nope" ) == -1 ); 
assert( cd["first"] == 1.0 ); 
```

The second assertion is the one that fails on 1.14.2 if the access is attempted, because the map
has gained a `{"nope", 0}` entry.
