#ifndef DFBA_CONSERVATION_LEDGER_H
#define DFBA_CONSERVATION_LEDGER_H

// Per-step conservation ledger for the extracellular substrates (R2.M1.7).
// Requires "./core/PhysiCell.h" to be included before this header.
//
// In this PhysiCell the cell sinks and sources are applied by phenotype.secretion.advance(), which
// update_all_cells() calls before the intracellular update; the diffusion solver
// (simulate_diffusion_decay) does not apply them. The ledger therefore has two windows:
//
//   cells window:  begin_cells(dt) before update_all_cells(), end_cells() after it.
//     dF = F(after) - F(before);  E = sum over active agents of net_export_rates[j]*dt, read on entry
//     (the rates that advance() applies);  residual = dF - E.
//   diffusion window:  begin_diffusion(dt) before simulate_diffusion_decay(), end_diffusion() after it.
//     dF = F(after) - F(before);  B = estimated diffusive inflow through Dirichlet voxels
//     (0 for a closed domain);  residual = dF - B.
//
// F_j = sum over voxels of density_j * voxel volume (mM*um^3). E uses the agents' net_export_rates,
// not the quantity BioFVM applies, so a stale cached rate or a lost update appears as a residual.
// "ext" is material added outside both windows (for example reinjecting oxygen): F at the start of a
// window minus F at the end of the previous one.
// Output: <folder>/conservation_ledger.csv, one row every `row_every` steps (cumulative sums).

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>
#include <omp.h>

namespace dfba_ledger
{

struct Ledger
{
	bool initialised = false;
	int nsub = 0;
	long long step = 0;
	int row_every = 100;
	FILE* fp = NULL;
	double dt_c = 0.0, dt_d = 0.0, t_last = 0.0;

	std::vector<double> F_c0, F_d0, F_last, E_dt, gross_dt, B0;
	std::vector<long double> c_dF, c_E, c_gross, c_res, d_dF, d_B, d_res, ext;
	std::vector<double> c_max, d_max;
	std::vector< std::pair<int,int> > bpairs;   // (Dirichlet voxel, interior neighbour)
	std::vector<int> bsub;

	~Ledger() { if( fp ) { write_row(); fclose( fp ); } }

	void total_content( std::vector<double>& F )
	{
		F.assign( nsub , 0.0 );
		long long nv = microenvironment.number_of_voxels();
		for( int j = 0 ; j < nsub ; j++ )
		{
			double s = 0.0;
			#pragma omp parallel for reduction(+:s)
			for( long long i = 0 ; i < nv ; i++ )
			{ s += microenvironment.density_vector(i)[j] * microenvironment.mesh.voxels[i].volume; }
			F[j] = s;
		}
	}

	// diffusive inflow rate (mM*um^3 per min) through the Dirichlet faces, from the current densities:
	// face flux between neighbouring voxels of volume V and spacing dx is D*(V/dx^2)*(rho_i - rho_k)
	void boundary_inflow( std::vector<double>& B )
	{
		B.assign( nsub , 0.0 );
		double dx = microenvironment.mesh.dx;
		for( size_t p = 0 ; p < bpairs.size() ; p++ )
		{
			int i = bpairs[p].first, k = bpairs[p].second, j = bsub[p];
			B[j] += microenvironment.diffusion_coefficients[j] * microenvironment.mesh.voxels[i].volume / ( dx * dx )
			        * ( microenvironment.density_vector(i)[j] - microenvironment.density_vector(k)[j] );
		}
	}

	void init()
	{
		nsub = microenvironment.number_of_densities();
		long long nv = microenvironment.number_of_voxels();
		for( long long i = 0 ; i < nv ; i++ )
		{
			if( microenvironment.is_dirichlet_node( (int) i ) == false ) { continue; }
			for( size_t c = 0 ; c < microenvironment.mesh.connected_voxel_indices[i].size() ; c++ )
			{
				int k = microenvironment.mesh.connected_voxel_indices[i][c];
				if( microenvironment.is_dirichlet_node( k ) ) { continue; }
				for( int j = 0 ; j < nsub ; j++ )
				{
					if( microenvironment.get_substrate_dirichlet_activation( j , (int) i ) )
					{ bpairs.push_back( std::make_pair( (int) i , k ) ); bsub.push_back( j ); }
				}
			}
		}
		c_dF.assign( nsub , 0.0 ); c_E.assign( nsub , 0.0 ); c_gross.assign( nsub , 0.0 ); c_res.assign( nsub , 0.0 );
		d_dF.assign( nsub , 0.0 ); d_B.assign( nsub , 0.0 ); d_res.assign( nsub , 0.0 ); ext.assign( nsub , 0.0 );
		c_max.assign( nsub , 0.0 ); d_max.assign( nsub , 0.0 );
		std::string name = PhysiCell_settings.folder + "/conservation_ledger.csv";
		fp = fopen( name.c_str() , "w" );
		if( fp )
		{
			fprintf( fp , "step,t_min" );
			for( int j = 0 ; j < nsub ; j++ )
			{
				const char* n = microenvironment.density_names[j].c_str();
				fprintf( fp , ",%s_F,%s_cells_dF,%s_cells_E,%s_cells_gross,%s_cells_resid,%s_cells_max_step_resid,%s_diff_dF,%s_diff_B,%s_diff_resid,%s_diff_max_step_resid,%s_ext" ,
					n,n,n,n,n,n,n,n,n,n,n );
			}
			fprintf( fp , "\n" );
		}
		initialised = true;
	}

	void write_row()
	{
		if( !fp || F_last.size() != (size_t) nsub ) { return; }
		fprintf( fp , "%lld,%.6f" , step , t_last );
		for( int j = 0 ; j < nsub ; j++ )
		{
			fprintf( fp , ",%.17g,%.17Lg,%.17Lg,%.17Lg,%.17Lg,%.17g,%.17Lg,%.17Lg,%.17Lg,%.17g,%.17Lg" ,
				F_last[j] , c_dF[j] , c_E[j] , c_gross[j] , c_res[j] , c_max[j] , d_dF[j] , d_B[j] , d_res[j] , d_max[j] , ext[j] );
		}
		fprintf( fp , "\n" ); fflush( fp );
	}

	void note_gap( const std::vector<double>& F_now )
	{
		if( F_last.size() == (size_t) nsub )
		{ for( int j = 0 ; j < nsub ; j++ ) { ext[j] += F_now[j] - F_last[j]; } }
	}
};

inline Ledger& get() { static Ledger L; return L; }

// ---- cells window (wrap update_all_cells) ----
inline void begin_cells( double dt )
{
	Ledger& L = get();
	if( !L.initialised ) { L.init(); }
	L.dt_c = dt;
	L.total_content( L.F_c0 );
	L.note_gap( L.F_c0 );
	L.E_dt.assign( L.nsub , 0.0 ); L.gross_dt.assign( L.nsub , 0.0 );
	long long na = all_basic_agents.size();
	for( int j = 0 ; j < L.nsub ; j++ )
	{
		double e = 0.0, g = 0.0;
		#pragma omp parallel for reduction(+:e,g)
		for( long long a = 0 ; a < na ; a++ )
		{
			Basic_Agent* p = all_basic_agents[a];
			if( p->is_active == false ) { continue; }
			// a newborn's net_export_rates pointer is bound to its phenotype only in its first advance(), so read the phenotype
			PhysiCell::Cell* pc = dynamic_cast<PhysiCell::Cell*>( p );
			double x = ( pc ? pc->phenotype.secretion.net_export_rates[j] : (*(p->net_export_rates))[j] ) * dt;
			e += x; g += std::fabs( x );
		}
		L.E_dt[j] = e; L.gross_dt[j] = g;
	}
}

inline void end_cells()
{
	Ledger& L = get();
	std::vector<double> F;
	L.total_content( F );
	for( int j = 0 ; j < L.nsub ; j++ )
	{
		double dF = F[j] - L.F_c0[j], r = dF - L.E_dt[j];
		L.c_dF[j] += dF; L.c_E[j] += L.E_dt[j]; L.c_gross[j] += L.gross_dt[j]; L.c_res[j] += r;
		if( std::fabs( r ) > L.c_max[j] ) { L.c_max[j] = std::fabs( r ); }
	}
	L.F_last = F;
}

// ---- diffusion window (wrap simulate_diffusion_decay) ----
inline void begin_diffusion( double dt )
{
	Ledger& L = get();
	if( !L.initialised ) { L.init(); }
	L.dt_d = dt;
	L.total_content( L.F_d0 );
	L.note_gap( L.F_d0 );
	L.boundary_inflow( L.B0 );
}

inline void end_diffusion( double t_after )
{
	Ledger& L = get();
	std::vector<double> F, B1;
	L.total_content( F );
	L.boundary_inflow( B1 );
	for( int j = 0 ; j < L.nsub ; j++ )
	{
		double dF = F[j] - L.F_d0[j], B = 0.5 * ( L.B0[j] + B1[j] ) * L.dt_d, r = dF - B;
		L.d_dF[j] += dF; L.d_B[j] += B; L.d_res[j] += r;
		if( std::fabs( r ) > L.d_max[j] ) { L.d_max[j] = std::fabs( r ); }
	}
	L.F_last = F;
	L.t_last = t_after;
	L.step++;
	if( L.step % L.row_every == 0 ) { L.write_row(); }
}

} // namespace dfba_ledger

#endif
