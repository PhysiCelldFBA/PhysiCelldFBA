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
//#include "../addons/dFBA/src/dfba_intracellular.h"


void create_cell_types(void)
{
	SeedRandom(parameters.ints("random_seed"));

	initialize_default_cell_definition();

	/*  This parses the cell definitions in the XML config file.  */
	initialize_cell_definitions_from_pugixml();

	//  This sets the pre and post intracellular update functions
	cell_defaults.functions.pre_update_intracellular =  NULL;
	cell_defaults.functions.post_update_intracellular = NULL;
	cell_defaults.functions.update_phenotype = NULL; 
	cell_defaults.functions.volume_update_function = NULL;

	build_cell_definitions_maps();
	
	setup_signal_behavior_dictionaries();

	Cell_Definition* c_beijerinckii = find_cell_definition( "c_beijerinckii");
	//  This sets the pre and post intracellular update functions for C. beijerinckii
	c_beijerinckii->functions.pre_update_intracellular =  NULL;
	c_beijerinckii->functions.post_update_intracellular = post_update_intracellular_c_beijerinckii;
	c_beijerinckii->functions.update_phenotype = NULL; 
	c_beijerinckii->functions.volume_update_function = NULL;

	Cell_Definition* m_barkeri = find_cell_definition( "m_barkeri");
	//  This sets the pre and post intracellular update functions for M. barkeri
	m_barkeri->functions.pre_update_intracellular =  NULL;
	m_barkeri->functions.post_update_intracellular = post_update_intracellular_m_barkeri;
	m_barkeri->functions.update_phenotype = NULL; 
	m_barkeri->functions.volume_update_function = NULL;

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
	// std::cout << std::endl; 
	
	// load cells from your CSV file
	load_cells_from_pugixml();
	
	return; 
}

// Cell-type specific intracellular update function for C. beijerinckii
void post_update_intracellular_c_beijerinckii(PhysiCell::Cell* pCell, PhysiCell::Phenotype& phenotype, double dt) {
	// Set growth rate from intracellular model
	pCell->custom_data["growth_rate"] = pCell->phenotype.intracellular->get_growth_rate();
	
	// Debug: Print growth rate and biomass flux for early timepoints
	if (PhysiCell_globals.current_time < 2.0) {
		double biomass_flux = 0.0;
		try {
			biomass_flux = pCell->phenotype.intracellular->get_flux_value("R_biomass");
		} catch (const std::exception& e) {
			std::cout << "Error getting biomass flux for c_beijerinckii: " << e.what() << std::endl;
		}
		
		// std::cout << "Time: " << PhysiCell_globals.current_time 
				//   << ", C. beijerinckii ID: " << pCell->ID 
				//   << ", Growth rate: " << pCell->custom_data["growth_rate"] 
				//   << ", Biomass flux: " << biomass_flux << std::endl;
	}
	
	// Set metabolic fluxes specific to C. beijerinckii
	try {
		pCell->custom_data["glucose_flux"] = pCell->phenotype.intracellular->get_flux_value("R_EX_glc_D_e");
		pCell->custom_data["hydrogen_flux"] = pCell->phenotype.intracellular->get_flux_value("R_EX_h2_e");
		
		// Debug fluxes for early timepoints (commented out for performance)
		// if (PhysiCell_globals.current_time < 5.0) {
		// 	std::cout << "C.beijerinckii ID " << pCell->ID 
		// 			  << ": Glucose=" << pCell->custom_data["glucose_flux"] 
		// 			  << ", H2=" << pCell->custom_data["hydrogen_flux"] << std::endl;
		// }
	} catch (const std::exception& e) {
		std::cout << "Error accessing flux values for c_beijerinckii: " << e.what() << std::endl;
	}
	
	return;
}

// Cell-type specific intracellular update function for M. barkeri
void post_update_intracellular_m_barkeri(PhysiCell::Cell* pCell, PhysiCell::Phenotype& phenotype, double dt) {
	// Set growth rate from intracellular model
	pCell->custom_data["growth_rate"] = pCell->phenotype.intracellular->get_growth_rate();
	
	// Debug: Print growth rate and biomass flux for early timepoints
	if (PhysiCell_globals.current_time < 2.0) {
		double biomass_flux = 0.0;
		try {
			biomass_flux = pCell->phenotype.intracellular->get_flux_value("R_BIOMASS_Mb_30");
		} catch (const std::exception& e) {
			std::cout << "Error getting biomass flux for m_barkeri: " << e.what() << std::endl;
		}
		
		// std::cout << "Time: " << PhysiCell_globals.current_time 
		// 		  << ", M. barkeri ID: " << pCell->ID 
		// 		  << ", Growth rate: " << pCell->custom_data["growth_rate"] 
		// 		  << ", Biomass flux: " << biomass_flux << std::endl;
	}
	
	// Set metabolic fluxes specific to M. barkeri
	try {
		pCell->custom_data["h2_flux"] = pCell->phenotype.intracellular->get_flux_value("R_EX_h2_e");
		pCell->custom_data["co2_flux"] = pCell->phenotype.intracellular->get_flux_value("R_EX_co2_e");
		pCell->custom_data["methane_flux"] = pCell->phenotype.intracellular->get_flux_value("R_EX_ch4_e");
		
		// Debug fluxes for early timepoints (commented out for performance)
		// if (PhysiCell_globals.current_time < 10.0) {
		// 	std::cout << "M.barkeri ID " << pCell->ID 
		// 			  << ": H2=" << pCell->custom_data["h2_flux"]
		// 			  << ", CO2=" << pCell->custom_data["co2_flux"]
		// 			  << ", CH4=" << pCell->custom_data["methane_flux"] << std::endl;
		// }
	} catch (const std::exception& e) {
		std::cout << "Error accessing flux values for m_barkeri: " << e.what() << std::endl;
	}
	
	return;
}

// Generic post_update_intracellular function that dispatches to cell-type specific functions
void post_update_intracellular(PhysiCell::Cell* pCell, PhysiCell::Phenotype& phenotype, double dt) {
	// Dispatch to the appropriate cell-type specific function
	if (pCell->type_name == "c_beijerinckii") {
		post_update_intracellular_c_beijerinckii(pCell, phenotype, dt);
	} else if (pCell->type_name == "m_barkeri") {
		post_update_intracellular_m_barkeri(pCell, phenotype, dt);
	} else {
		// Default behavior for unknown cell types
		std::cout << "Warning: Unknown cell type '" << pCell->type_name 
				  << "' in post_update_intracellular" << std::endl;
		pCell->custom_data["growth_rate"] = pCell->phenotype.intracellular->get_growth_rate();
	}
	
	return;
}


// NOT USED
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
        
	double flux_value =  0.0;
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
