#include "dfba_intracellular.h"
#include <sstream>
#include <iostream>

namespace PhysiCelldFBA {

dFBAIntracellular::dFBAIntracellular() : Intracellular()
{
    this->intracellular_type = "dfba";
}

dFBAIntracellular::dFBAIntracellular(pugi::xml_node& node)
{
    intracellular_type = "dfba";
	this->initialize_intracellular_from_pugixml(node);
}

dFBAIntracellular::dFBAIntracellular(dFBAIntracellular* copy) 
{
    this->intracellular_type = copy->intracellular_type;
    this->sbml_filename = copy->sbml_filename;
    this->cell_density = copy->cell_density;
    this->max_growth_rate = copy->max_growth_rate;
    this->current_growth_rate = copy->current_growth_rate;
    this->next_dfba_run = copy->next_dfba_run;
    //this->model = new dFBAModel();
    this->model.initModel(this->sbml_filename.c_str());
}


int dFBAIntracellular::parse_transport_model(pugi::xml_node& node)
{
    ExchangeFluxData exchange_flux;
    KineticParam Km;
    KineticParam Vmax;
    
    pugi::xml_node node_exchange = node.child( "exchange" );
	int num_exchanges = 0;
    while( node_exchange )
	{
		string density_name = node_exchange.attribute( "substrate" ).value(); 
        int density_index = microenvironment.find_density_index( density_name ); 
        std::cout << "Parsing " << density_name << std::endl;
        std::string actual_name = microenvironment.density_names[ density_index ]; 
			
        // error check 
        if( std::strcmp( density_name.c_str() , actual_name.c_str() ) != 0 )
        {
            std::cout << "Error: attempted to set secretion/uptake/export for \"" 
                << density_name << "\", which was not found in the microenvironment." << std::endl 
            << "       Please double-check your substrate name in the config file." << std::endl << std::endl; 
            exit(-1); 
        }
        
        pugi::xml_node node_fba_flux = node_exchange.child( "fba_flux" ); 
		if( node_fba_flux )
		{  
            exchange_flux.fba_flux_id = PhysiCell::xml_get_my_string_value(node_fba_flux);
        }
        else {
            std::cout << "Error: attempted get fba_flux node for "; 
            std::cout << exchange_flux.density_name << "\", but not found." << std::endl;
            std::cout << "Please double-check your exchange nodes in the config file." << std::endl;
            std::cout << std::endl; 
            exit(-1); 
        }

        pugi::xml_node node_Km = node_exchange.child( "Km" ); 
		if( node_Km )
		{
            Km.name = "Km";
            Km.untis = node_Km.attribute("units").value();
            Km.value = PhysiCell::xml_get_my_double_value(node_Km);
            
        }
        else {
            std::cout << "Error: attempted get Km node for "; 
            std::cout << density_name << "\", but not found." << std::endl;
            std::cout << "Please double-check your exchange nodes in the config file." << std::endl;
            std::cout << std::endl; 
            exit(-1); 
        }

        pugi::xml_node node_Vmax = node_exchange.child( "Vmax" ); 
		if( node_Vmax )
		{
            Vmax.name = "Vmax";
            Vmax.untis = node_Vmax.attribute("units").value();
            Vmax.value = PhysiCell::xml_get_my_double_value(node_Vmax);
            
        }
        else {
            std::cout << "Error: attempted get Vmax node for "; 
            std::cout << density_name << "\", but not found." << std::endl;
            std::cout << "Please double-check your exchange nodes in the config file." << std::endl;
            std::cout << std::endl; 
            exit(-1); 
        }
		
        exchange_flux.density_name = density_name;
        exchange_flux.density_index = density_index;
        exchange_flux.Km = Km;
        exchange_flux.Vmax = Vmax;

        this->substrate_exchanges[density_name] = exchange_flux;
        num_exchanges++;
		node_exchange = node_exchange.next_sibling( "exchange" ); 
	}
    return num_exchanges;
}

void dFBAIntracellular::parse_growth_model(pugi::xml_node& parent)
{
    pugi::xml_node node = parent.child( "cell_density" );
	if ( node )
	{ 
        this->cell_density = PhysiCell::xml_get_my_double_value(node);
    }
    else
    {
        std::cout << "Error: attempted to read sbml_filename path but not found." << std::endl;
        std::cout << "Please double-check your exchange nodes in the XML setting." << std::endl;
        std::cout << std::endl; 
        exit(-1); 
    }
}

void dFBAIntracellular::initialize_intracellular_from_pugixml(pugi::xml_node& node)
{

    // Getting sbml file name for reading the model
    pugi::xml_node node_sbml = node.child( "sbml_filename" );
	if ( node_sbml )
	{ 
        this->sbml_filename = PhysiCell::xml_get_my_string_value (node_sbml);
    }
    else
    {
        std::cout << "Error: attempted to read sbml_filename path but not found." << std::endl;
        std::cout << "Please double-check your exchange nodes in the XML setting." << std::endl;
        std::cout << std::endl; 
        exit(-1); 
    }
	

    // parsing the transport model
    pugi::xml_node node_transport_model = node.child( "transport_model" );
    if ( node_transport_model )
	{ 
        int num_exchanges = parse_transport_model(node_transport_model);
        if (num_exchanges == 0){
            std::cout << "Error, dFBA model must have at least one exchange flux." << std::endl;
            std::cout << "Please double-check your exchange nodes in the XML setting file." << std::endl;
            std::cout << std::endl; 
            exit(-1);
        }
    }
    else
    {
        std::cout << "Error: attempted to parse transport_model but not found." << std::endl;
        std::cout << "Please double-check intracellular model nodes in the XML setting file." << std::endl;
        std::cout << std::endl; 
        exit(-1); 
    }


    // parsing the transport model
    pugi::xml_node node_growth_model = node.child( "growth_model" );
    if ( node_growth_model )
	{ 
        parse_growth_model(node_growth_model);
    }
    else
    {
        std::cout << "Error: attempted to parse growth_model but not found." << std::endl;
        std::cout << "Please double-check intracellular model nodes in the XML setting file." << std::endl;
        std::cout << std::endl; 
        exit(-1); 
    }
}



void dFBAIntracellular::start()
{

    std::cout << "Loaing SBML model from: " << this->sbml_filename << std::endl;
    dFBAModel metnet;
    metnet.readSBMLModel(this->sbml_filename.c_str());

    int col_idx = 0;
    for(dFBAReaction* rxn: metnet.getListOfReactions()){
        this->reactionsIndexer[rxn->getId()] = col_idx;
        col_idx++;
    }
    
    int row_index = 0;
    for( dFBAMetabolite* met: metnet.getListOfMetabolites() ){
        this->metaboliteIndexer[met->getId()] = row_index;
        row_index++;
    }

    this->handler = new CoinMessageHandler(nullptr);
    this->handler->setLogLevel(0);
    this->problem.passInMessageHandler(handler);

    int n_rows = metnet.getNumMetabolites();
    int n_cols = metnet.getNumReactions();

    CoinPackedMatrix matrix;
    matrix.setDimensions(n_rows, 0);

    double* row_lb = new double[n_rows]; //the row lower bounds
    double* row_ub = new double[n_rows]; //the row upper bounds
    double* col_lb = new double[n_cols]; //the column lower bounds
    double* col_ub = new double[n_cols]; //the column upper bounds
    double* objective = new double[n_cols]; //the objective coefficients

    for(int i=0; i< n_rows; i++)
    {
        row_lb[i] = 0;
        row_ub[i] = 0;
    }

    for(dFBAReaction* rxn: metnet.getListOfReactions())
    {
        int col_idx = this->reactionsIndexer[rxn->getId()];
       
        col_lb[col_idx] = rxn->getLowerBound();
        col_ub[col_idx] = rxn->getUpperBound();
        objective[col_idx] = rxn->getObjectiveCoefficient();

        const std::map<const dFBAMetabolite*, double> metabolites = rxn->getMetabolites();
        CoinPackedVector col;
        for(auto it=metabolites.begin(); it!=metabolites.end(); it++)
        {
            const dFBAMetabolite* met = it->first;
            double stoich_coeff = it->second;
            int row_idx = this->metaboliteIndexer[met->getId()];
            col.insert(row_idx, stoich_coeff);
        }
        matrix.appendCol(col);
    }

    this->problem.loadProblem(matrix, col_lb, col_ub, objective, row_lb, row_ub);
    this->problem.setOptimizationDirection(-1);

    delete col_lb;
    delete col_ub;
    delete row_lb;
    delete row_ub;
    delete objective;

    this->is_initialized = true;
}


void update_dfba_inputs( PhysiCell::Cell* pCell, PhysiCell::Phenotype& phenotype, double dt ){

void dFBAIntracellular::update(PhysiCell::Cell* pCell, PhysiCell::Phenotype& phenotype, double dt)
{
    bool debug = false;

    /*
        Steps for the update
        1- Update exchange fluxes lower bound using current concentrations
        2- run FBA and check growth threshold
        3- update the cell volumne using the growth rate from FBA
        4- rescale exchange fluxes from the dfba model and use them to update the net_export_rates
        5- remove the internalized substrates if needed

        mM /gDW cell / h <-- fluxes uints
        Biomass: 

    */

    // Setp 1 - Update exchange fluxes lower
    map<std::string, exchange_data>::iterator it;
    for(it = this->substrate_exchanges.begin(); it != this->substrate_exchanges.end(); it++)
    {
        std::string substrate_name = it->first;
        ExchangeFluxData ex_strut = it->second;

        // geting the amount of substrate
        double substrate_conc = density_vector[ex_strut.density_index];;
        // scaling Vmax ased on cell volume
        double Vmax = ex_strut.Vmax.value * Vmax_scale;
        double Km = ex_strut.Km.value;
        
        // geting the amount of substrate
        double substrate_conc = pCell->nearest_density_vector()[ex_strut.density_index];

        // PROBLEM TO SOLVE
        // having substrate_conc in the proper units that match FBA

        // useing irreversible Michaelis Menten kinetics to estimate the flux bound
        double flux_bound = (Vmax * substrate_conc) / (Km + substrate_conc); // should be calculated from density
        // Change sign to use as lower bound of the exchange flux
        double exchange_flux_lb = -1 * max_rate;
        // Updateing the lower bound of the corresponding exchange flux
        this->model.setReactionLowerBound(ex_strut.fba_flux_id, flux_bound);    
        
        if ( debug ) {
            std::cout << " - [" << substrate_name << "] = " << substrate_conc;
            std::cout << " ==> " << ex_strut.fba_flux_id << " = " << flux_bound << std::endl;
        }
    }

    // STEP 2 - run FBA
    this->model.runFBA();

    // STEP 3 - 
    float growth_rate = this->model.getObjectiveValue();
    
    // V(t+1) = V(t) + V(t) * u * dt = V(t) * (1 + u * dt)
    // growth_rate 1/h
    float volume_increase_ratio = 1 + (growth_rate / 60) * dt
    phenotype.volume.multiply_by_ratio( volume_increase_ratio );
    
    
    // STEPS 4-5
    map<std::string, exchange_data>::iterator it;
    for(it = this->substrate_exchanges.begin(); it != this->substrate_exchanges.end(); it++)
    {
        // Retrive the exchange flux and its corresponding density
        std::string substrate_name = it->first;
        ExchangeFluxData ex_strut = it->second;
        
        int density_index = ex_strut.density_index;
        std::string fba_flux_id = ex_strut.fba_flux_id;
        
        FBA_reaction* exchange = this->model.getReaction(fba_flux_id);
        // flux units: mmol / gDW cell / h
        double flux = exchange->getFluxValue();
        flux *= -1;
        
        // voxel volumne = 8000 um³ = 0.8e-11 L
        // FLUXES Units: mmol / gDW cell / hr
        // µm3 = 1 fL = 1e−15 L
        // MULTIPLYER (hr) --> * (dt/60)

        // HeLa cell mass 2.3 ng
        // cell_density --> 1.04 g/ml = 1.04 ug/nL = 1.04 ng / pL = 0.00104 ng / fL
        // cell volume fL (um³)
        
        // cell.mass.total (ng) = cell.volume (um³) * cell.density (ng / um³)
        // ~2.6 (ng) = ~2500 (um³) * ~0.00104 ng / (um³)
        
        // gDW cell (cell.volume.total * cell_density) = mass.total
        // cell.mass.solid = cell.mass.total * (1-fluid_frac)

        // MULTIPLYER (gDW cell) --> * cell.mass.solid (ng DW cell)
        
        
        // mmol / per_gDW_per_hr 
        
        // 
        // mM: mmol / L
        // how to rescale FBA exchanges into net_export_rates
        // ml = 1e+12 um³
        // g = 1e+12 pg
        // cell density (from HeLa) g / ml = pg / um³
        
        float cell_density = 1.04;
        float solid_fraction = (1 - phenotype.volume.fluid_fraction) ;
        // dry wight units: pg
        float dry_weight = (phenotype.volume.total * cell_density) * solid_fraction
        float scaling = 1;
        
        float net_export_rate = flux * dry_weight * dt
        
        flux *= scaling;
        phenotype.secretion.net_export_rates[density_index] = flux;

        // STEP 4 - setting internalized total substrate to 0
        phenotype.molecular.internalized_total_substrates[density_index] = 0;

        
    }
    


}

int dFBAIntracellular::update_phenotype_parameters(PhysiCell::Phenotype& phenotype)
{


    return 0;
}

void dFBAIntracellular::save_dFBA(std::string path, std::string index) 
{
	
}

}