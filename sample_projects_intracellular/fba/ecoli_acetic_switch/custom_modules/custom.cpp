/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2018, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "custom.h"
#include <set>   // for the one-time warning in set_custom_data
//#include "../addons/dFBA/src/dfba_intracellular.h"


// ---------------------------------------------------------------------------
// Safe custom_data write.
//
// PhysiCell's Custom_Cell_Data::operator[](std::string) resolves a name via
// std::map::operator[], which INSERTS a missing key with value 0 -- so writing
// to a name that is not declared in the XML silently lands in variables[0] and
// corrupts whatever variable happens to be first. Use this instead: an
// undeclared name is reported once and skipped rather than clobbering data.
//
// The report is guarded because this runs on the OpenMP-parallel update path.
// ---------------------------------------------------------------------------
static void set_custom_data(PhysiCell::Cell* pCell, const std::string& name, double value)
{
	const int index = pCell->custom_data.find_variable_index(name);
	if (index < 0)
	{
		static std::set<std::string> reported;
		#pragma omp critical(dfba_custom_data_warn)
		{
			if (reported.insert(name).second)
			{
				std::cout << "WARNING: custom_data variable \"" << name
				          << "\" is written from custom.cpp but not declared in the"
				          << " XML cell definition; the write is being skipped."
				          << std::endl;
			}
		}
		return;
	}
	pCell->custom_data.variables[index].value = value;
}

void create_cell_types(void)
{
	SeedRandom(parameters.ints("random_seed"));

	initialize_default_cell_definition();

	/*  This parses the cell definitions in the XML config file.  */
	initialize_cell_definitions_from_pugixml();

	//  This sets the pre and post intracellular update functions
	cell_defaults.functions.pre_update_intracellular =  NULL;
	cell_defaults.functions.post_update_intracellular = post_update_intracellular;
	cell_defaults.functions.update_phenotype = NULL; 
	cell_defaults.functions.volume_update_function = NULL;

	build_cell_definitions_maps();
	
	setup_signal_behavior_dictionaries();

	Cell_Definition* ecoli = find_cell_definition( "ecoli");
	//  This sets the pre and post intracellular update functions
	ecoli->functions.pre_update_intracellular =  NULL;
	ecoli->functions.post_update_intracellular = post_update_intracellular;
	ecoli->functions.update_phenotype = NULL; 
	ecoli->functions.volume_update_function = NULL;

	display_cell_definitions(std::cout);

	return;
}




void setup_microenvironment(void)
{
	initialize_microenvironment();
	return;
}

void setup_tissue(void)
{

	double Xmin = microenvironment.mesh.bounding_box[0]; 
	double Ymin = microenvironment.mesh.bounding_box[1]; 
	double Zmin = microenvironment.mesh.bounding_box[2]; 

	double Xmax = microenvironment.mesh.bounding_box[3]; 
	double Ymax = microenvironment.mesh.bounding_box[4]; 
	double Zmax = microenvironment.mesh.bounding_box[5]; 
	
	if( default_microenvironment_options.simulate_2D == true )
	{
		Zmin = 0.0; 
		Zmax = 0.0; 
	}
	
	double Xrange = Xmax - Xmin; 
	double Yrange = Ymax - Ymin; 
	double Zrange = Zmax - Zmin; 
	
	// create some of each type of cell 
	
	Cell* pC;
	
	for( int k=0; k < cell_definitions_by_index.size() ; k++ )
	{
		Cell_Definition* pCD = cell_definitions_by_index[k]; 
		std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl; 
		for( int n = 0 ; n < parameters.ints("number_of_cells") ; n++ )
		{
			std::vector<double> position = {0,0,0}; 
			position[0] = Xmin + UniformRandom()*Xrange; 
			position[1] = Ymin + UniformRandom()*Yrange; 
			position[2] = Zmin + UniformRandom()*Zrange; 
			
			pC = create_cell( *pCD ); 
			pC->assign_position( position );
		}
	}
	std::cout << std::endl; 
	
	// load cells from your CSV file
	load_cells_from_pugixml();
	
	return; 
}

void post_update_intracellular(PhysiCell::Cell* pCell, PhysiCell::Phenotype& phenotype, double dt ){

	
	PhysiCelldFBA::dFBAIntracellular* dfba_model =
		static_cast<PhysiCelldFBA::dFBAIntracellular*>(phenotype.intracellular);

	// The project's established names, kept so existing analysis scripts
	// (scripts/dfba_analysis.py and the large-scale notebook) keep working.
	set_custom_data(pCell, "growth_rate", dfba_model->get_growth_rate());
	set_custom_data(pCell, "oxygen_flux",  dfba_model->get_flux_value("R_EX_o2_e"));
	set_custom_data(pCell, "glucose_flux", dfba_model->get_flux_value("R_EX_glc__D_e"));
	set_custom_data(pCell, "acetate_flux", dfba_model->get_flux_value("R_EX_ac_e"));
	set_custom_data(pCell, "co2_flux",     dfba_model->get_flux_value("R_EX_co2_e"));

	// Keyed by reaction id, for the FVA degeneracy diagnostic: the realised flux
	// and the Michaelis-Menten lower bound that constrained it. Together they
	// let the cell's LP be replayed offline. Only the lower bound is recorded:
	// update_dfba_inputs never touches an exchange's upper bound.
	for (const auto& exchange : dfba_model->substrate_exchanges)
	{
		const PhysiCelldFBA::ExchangeFluxData& ex = exchange.second;
		set_custom_data(pCell, ex.fba_flux_id,
			dfba_model->get_flux_value(ex.fba_flux_id));
		set_custom_data(pCell, ex.fba_flux_id + "_lb",
			dfba_model->sbml_model.getReactionLowerBound(ex.fba_flux_id));
	}

	return;
}

void reintroduce_nutrients_function () 
{
	if (PhysiCell::parameters.bools.find_index("nutrient_reintroduction") != -1) 
	{
		int nutrient_index = BioFVM::microenvironment.find_density_index(PhysiCell::parameters.strings("reintroduced_nutrient"));

		if (PhysiCell::parameters.bools("nutrient_reintroduction")){
			// Activate nutrient boundary at specified time
			if (
				(PhysiCell::PhysiCell_globals.current_time >= PhysiCell::parameters.doubles("reintroduction_start_time"))
				&& (PhysiCell::PhysiCell_globals.current_time < (PhysiCell::parameters.doubles("reintroduction_start_time") + PhysiCell::parameters.doubles("reintroduction_duration")))
				&& !BioFVM::microenvironment.get_substrate_dirichlet_activation(nutrient_index)
			)
			{
				std::cout << PhysiCell::parameters.strings("reintroduced_nutrient") << " boundary activated at t=" << PhysiCell::PhysiCell_globals.current_time << std::endl;
				BioFVM::microenvironment.set_substrate_dirichlet_activation(nutrient_index, true);	
				std::cout << "Boundary condition set at: " << BioFVM::microenvironment.get_substrate_dirichlet_value(nutrient_index, 0) << std::endl;
			}
			else if (PhysiCell::PhysiCell_globals.current_time < (PhysiCell::parameters.doubles("reintroduction_start_time")))
			{
				BioFVM::microenvironment.set_substrate_dirichlet_activation(nutrient_index, false);
			}

			// Deactivate nutrient boundary after the duration
			if (
				(PhysiCell::PhysiCell_globals.current_time >= (PhysiCell::parameters.doubles("reintroduction_start_time") + PhysiCell::parameters.doubles("reintroduction_duration")))
				&& BioFVM::microenvironment.get_substrate_dirichlet_activation(nutrient_index)
			)
			{
				std::cout << PhysiCell::parameters.strings("reintroduced_nutrient") << " boundary deactivated at t=" << PhysiCell::PhysiCell_globals.current_time << std::endl;
				BioFVM::microenvironment.set_substrate_dirichlet_activation(nutrient_index, false);	
			}
			
		} else if ( BioFVM::microenvironment.get_substrate_dirichlet_activation(nutrient_index) ){
			std::cout << PhysiCell::parameters.strings("reintroduced_nutrient") << " boundary forced deactivation at t=" << PhysiCell::PhysiCell_globals.current_time << std::endl;
			BioFVM::microenvironment.set_substrate_dirichlet_activation(nutrient_index, false);	
		}
	}
}

void inject_density(int density_index, double concentration)
{
	// Inject given concentration on the extremities only
	#pragma omp parallel for
	for (int n = 0; n < microenvironment.number_of_voxels(); n++)
	{
		auto current_voxel = microenvironment.voxels(n);
		microenvironment.density_vector(n)[density_index] = concentration;
	}
}



std::vector<std::vector<double>> create_cell_sphere_positions(double cell_radius, double sphere_radius)
{
	std::vector<std::vector<double>> cells;
	int xc=0,yc=0,zc=0;
	double x_spacing= cell_radius*sqrt(3);
	double y_spacing= cell_radius*2;
	double z_spacing= cell_radius*sqrt(3);
	
	std::vector<double> tempPoint(3,0.0);
	// std::vector<double> cylinder_center(3,0.0);
	
	for(double z=-sphere_radius;z<sphere_radius;z+=z_spacing, zc++)
	{
		for(double x=-sphere_radius;x<sphere_radius;x+=x_spacing, xc++)
		{
			for(double y=-sphere_radius;y<sphere_radius;y+=y_spacing, yc++)
			{
				tempPoint[0]=x + (zc%2) * 0.5 * cell_radius;
				tempPoint[1]=y + (xc%2) * cell_radius;
				tempPoint[2]=z;
				
				if(sqrt(norm_squared(tempPoint))< sphere_radius)
				{ cells.push_back(tempPoint); }
			}
			
		}
	}
	return cells;
	
}


std::vector<std::vector<double>> create_cell_disc_positions(double cell_radius, double disc_radius)
{	 
	double cell_spacing = 0.95 * 2.0 * cell_radius; 
	
	double x = 0.0; 
	double y = 0.0; 
	double x_outer = 0.0;

	std::vector<std::vector<double>> positions;
	std::vector<double> tempPoint(3,0.0);
	
	int n = 0; 
	while( y < disc_radius )
	{
		x = 0.0; 
		if( n % 2 == 1 )
		{ x = 0.5 * cell_spacing; }
		x_outer = sqrt( disc_radius*disc_radius - y*y ); 
		
		while( x < x_outer )
		{
			tempPoint[0]= x; tempPoint[1]= y;	tempPoint[2]= 0.0;
			positions.push_back(tempPoint);			
			if( fabs( y ) > 0.01 )
			{
				tempPoint[0]= x; tempPoint[1]= -y;	tempPoint[2]= 0.0;
				positions.push_back(tempPoint);
			}
			if( fabs( x ) > 0.01 )
			{ 
				tempPoint[0]= -x; tempPoint[1]= y;	tempPoint[2]= 0.0;
				positions.push_back(tempPoint);
				if( fabs( y ) > 0.01 )
				{
					tempPoint[0]= -x; tempPoint[1]= -y;	tempPoint[2]= 0.0;
					positions.push_back(tempPoint);
				}
			}
			x += cell_spacing; 
		}		
		y += cell_spacing * sqrt(3.0)/2.0; 
		n++; 
	}
	return positions;
}


std::vector<std::string> my_coloring_function( Cell* pCell )
{

	std::vector<std::string> output(4, "red");

    std::string fba_flux_id = "R_EX_o2_e";
        
	double flux_value =  pCell->phenotype.intracellular->get_flux_value(fba_flux_id);

	if( abs(flux_value) >= 0.0 )
	{
		output[0] = "blue";
		output[2] = "blue";
		return output;
	}
	else
	{
		output[0] = "green";
		output[2] = "green";
	}

	return output;
}
