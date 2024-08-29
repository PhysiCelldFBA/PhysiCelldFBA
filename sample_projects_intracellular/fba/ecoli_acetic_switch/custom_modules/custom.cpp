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


void create_cell_types(void)
{
	SeedRandom(parameters.ints("random_seed"));

	initialize_default_cell_definition();

	/*  This parses the cell definitions in the XML config file.  */
	initialize_cell_definitions_from_pugixml();

	//  This sets the pre and post intracellular update functions
	cell_defaults.functions.pre_update_intracellular =  pre_update_intracellular;
	cell_defaults.functions.post_update_intracellular = post_update_intracellular;
	cell_defaults.functions.update_phenotype = NULL; 
	

	build_cell_definitions_maps();
	
	setup_signal_behavior_dictionaries();

	display_cell_definitions(std::cout);

	return;
}




void setup_microenvironment( void )
{
	initialize_microenvironment();
	return;
}


void setup_tissue( void )
{
	// place a bacterial colony at the center 
	double cell_radius = cell_defaults.phenotype.geometry.radius; 
	double cell_spacing = 0.95 * 2.0 * cell_radius; 
	double colony_radius = parameters.doubles("colony_radius");
	
	Cell* pCell = NULL; 
	
	double x = 0.0; 
	double x_outer = colony_radius; 
	double y = 0.0; 
	
	int n = 0; 
	while( y < colony_radius )
	{
		x = 0.0; 
		if( n % 2 == 1 )
		{ x = 0.5*cell_spacing; }
		x_outer = sqrt( colony_radius*colony_radius - y*y ); 
		
		while( x < x_outer )
		{
			pCell = create_cell(bacteria_cell);
			pCell->assign_position( x , y , 0.0 );
			
			if( fabs( y ) > 0.01 )
			{
				pCell = create_cell(bacteria_cell);
				pCell->assign_position( x , -y , 0.0 );
				pCell->fba_model = FBA::FBA_default_model;			
			}
			
			if( fabs( x ) > 0.01 )
			{ 
				pCell = create_cell(bacteria_cell);
				pCell->assign_position( -x , y , 0.0 );
				pCell->fba_model = FBA::FBA_default_model;			
				if( fabs( y ) > 0.01 )
				{
					pCell = create_cell(bacteria_cell);
					pCell->assign_position( -x , -y , 0.0 );
					pCell->fba_model = FBA::FBA_default_model;			
				}
			}
			x += cell_spacing; 
			
		}
		
		y += cell_spacing * sqrt(3.0)/2.0; 
		n++; 
	}
	
	return; 
}

void pre_update_intracellular(PhysiCell::Cell* pCell, PhysiCell::Phenotype& phenotype, double dt ){
	PhysiCelldFBA::update_dfba_inputs( pCell, phenotype, dt );
	PhysiCelldFBA::
	return;
}


void post_update_intracellular(PhysiCell::Cell* pCell, PhysiCell::Phenotype& phenotype, double dt ){
	PhysiCelldFBA::update_dfba_outputs( pCell, phenotype, dt );
	return;
}




void update_cell(PhysiCell::Cell* pCell, PhysiCell::Phenotype& phenotype, double dt ){

  phenotype.secretion.set_all_secretion_to_zero();
  phenotype.secretion.set_all_uptake_to_zero();

  std::string oxygen_name = "oxygen";
  std::string glucose_name = "glucose";
  std::string acetate_name = "acetate";
  
  std::string oxygen_flux_id = FBA::exchange_flux_density_map[oxygen_name];
  std::string glucose_flux_id = FBA::exchange_flux_density_map[glucose_name];
  std::string acetate_flux_id = FBA::exchange_flux_density_map[acetate_name];

  static int oxygen_idx = microenvironment.find_density_index( oxygen_name );
  static int glucose_idx = microenvironment.find_density_index( glucose_name );
  static int acetate_idx = microenvironment.find_density_index( acetate_name );

  double oxygen_density = pCell->nearest_density_vector()[oxygen_idx];
  double glucose_density = pCell->nearest_density_vector()[glucose_idx]; // dived by voxel size?
  double acetate_density = pCell->nearest_density_vector()[acetate_idx]; // dived by voxel size?
  
  
  double O2_Km = parameters.doubles("oxygen_Km");
  double O2_Vmax = parameters.doubles("oxygen_Vmax");
  double oxygen_flux_bound = -1 * ( O2_Vmax * oxygen_density) / (oxygen_density + O2_Km);
  pCell->fba_model.setReactionLowerBound(oxygen_flux_id, oxygen_flux_bound);

  double glc_Km = parameters.doubles("glucose_Km");
  double glc_Vmax = parameters.doubles("glucose_Vmax");
  double glucose_flux_bound = -1 * ( glc_Vmax * glucose_density) / (glucose_density + glc_Km);
  pCell->fba_model.setReactionLowerBound(glucose_flux_id, glucose_flux_bound);

  double lac_Km = parameters.doubles("acetate_Km");
  double lac_Vmax = parameters.doubles("acetate_Vmax");
  double acetate_flux_bound = -1 * ( lac_Vmax * acetate_density) / (acetate_density + lac_Km);
  pCell->fba_model.setReactionLowerBound(acetate_flux_id, acetate_flux_bound);

  std::cout << "Running FBA for cell: " << pCell->ID << std::endl;
  pCell->fba_model.runFBA();

  if ( pCell->fba_model.getSolutionStatus() )
  {
    FBA_reaction* reaction = pCell->fba_model.getReaction(oxygen_flux_id);
    double oxygen_flux = reaction->getFluxValue();
	std::cout << "Oxygen flux: " << oxygen_flux << std::endl;

	reaction = pCell->fba_model.getReaction(glucose_flux_id);
	double glucose_flux = reaction->getFluxValue();
	std::cout << "glucose flux: " << glucose_flux << std::endl;

	reaction = pCell->fba_model.getReaction(acetate_flux_id);
	double acetate_flux = reaction->getFluxValue();
	std::cout << "acetate flux: " << acetate_flux << std::endl;
		
	if ( oxygen_flux < 0)
		phenotype.secretion.uptake_rates[oxygen_idx] = abs(oxygen_flux / oxygen_density);

	if ( glucose_flux < 0)
		phenotype.secretion.uptake_rates[glucose_idx] = abs(glucose_flux / glucose_density);

	if ( acetate_flux < 0 )
		phenotype.secretion.uptake_rates[acetate_idx] = abs(acetate_flux / acetate_density);
			
	else if ( acetate_flux > 0 )
		phenotype.secretion.secretion_rates[acetate_idx] = abs(acetate_flux / acetate_density);
		
  }
  else
  {
	// Check energy production to see is the cell is able
	// to reach a threshold. Otherwise enter apoptosis

  }


    
  


}

void setup_default_metabolic_model( void )
{

  // This codes is to use iSIM model
  std::string sbml_fileame = parameters.strings("sbml_model");
  
  FBA::FBA_default_model.readSBMLModel(sbml_fileame.c_str());
  FBA::FBA_default_model.initLpModel();

  FBA::exchange_flux_density_map["glucose"] = "R_EX_glc__D_e";
  FBA::exchange_flux_density_map["acetate"] = "R_EX_ac_e";
  FBA::exchange_flux_density_map["oxygen"] = "R_EX_o2_e";
  FBA::exchange_flux_density_map["growth_rate"] = "R_BIOMASS_Ecoli_core_w_GAM";

}

void anuclear_volume_model (Cell* pCell, Phenotype& phenotype, double dt)
{
	phenotype.volume.fluid += dt * phenotype.volume.fluid_change_rate * 
		( phenotype.volume.target_fluid_fraction * phenotype.volume.total - phenotype.volume.fluid );
		
	// if the fluid volume is negative, set to zero
	if( phenotype.volume.fluid < 0.0 )
	{ phenotype.volume.fluid = 0.0; }
	

	phenotype.volume.cytoplasmic_solid += dt * phenotype.volume.cytoplasmic_biomass_change_rate * 
		( phenotype.volume.target_solid_cytoplasmic - phenotype.volume.cytoplasmic_solid );	
	
	if( phenotype.volume.cytoplasmic_solid < 0.0 )
	{ phenotype.volume.cytoplasmic_solid = 0.0; }
	
	phenotype.volume.solid = phenotype.volume.cytoplasmic_solid;
	phenotype.volume.cytoplasmic = phenotype.volume.cytoplasmic_solid + phenotype.volume.cytoplasmic_fluid; 
	
	phenotype.volume.total = phenotype.volume.cytoplasmic;
	
	
	phenotype.volume.fluid_fraction = phenotype.volume.fluid / 
		( 1e-16 + phenotype.volume.total ); 
   
	phenotype.geometry.update( pCell,phenotype,dt );

	return; 
}

void metabolic_cell_phenotype( Cell* pCell, Phenotype& phenotype, double dt )
{
	// if cell is dead, don't bother with future phenotype changes.
	if( phenotype.death.dead == true )
	{
		pCell->functions.update_phenotype = NULL;
		return;
	}

	// update the transition rate according to growth rate?
	static int cycle_start_index = live.find_phase_index( PhysiCell_constants::live );
	static int cycle_end_index = live.find_phase_index( PhysiCell_constants::live );

	//static int oncoprotein_i = pCell->custom_data.find_variable_index( "oncoprotein" );
	//phenotype.cycle.data.transition_rate( cycle_start_index ,cycle_end_index ) *= pCell->custom_data[oncoprotein_i] ;
	return;
}



std::vector<std::string> my_coloring_function( Cell* pCell )
{
	// start with flow cytometry coloring

	std::vector<std::string> output = false_cell_coloring_cytometry(pCell);
	output[0] = "red";
	output[1] = "red";
	output[2] = "red";

	if( pCell->phenotype.death.dead == false && pCell->type == 1 )
	{
		 output[0] = "black";
		 output[2] = "black";
	}

	return output;
}
