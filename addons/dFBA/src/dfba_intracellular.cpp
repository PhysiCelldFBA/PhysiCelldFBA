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
        std::string actual_name = microenvironment.density_names[ density_index ]; 
			
        // error check 
        if( std::strcmp( density_name.c_str() , actual_name.c_str() ) != 0 )
        {
            std::cout << "Error: attempted to set secretion/uptake/export for \""; 
            std::cout << density_name << "\", which was not found in the microenvironment." << std::endl;
            std::cout << "Please double-check your substrate name in the config file." << std::endl;
            std::cout << std::endl;
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

    // HeLa cell mass 2.3 ng
    // cell volume fL (um³)
    // mM: mmol / L = 10⁻³ mol / 1e¹⁵ um³ = 10⁻¹⁵mol / um³ = amol /  um³
    // cell_density --> 1.04 g/ml = 1.04 ug/nL = 1.04 ng / pL = 1.04 pg / fL = pg / um³
    // cell.mass.total (ng) = cell.volume (um³) * cell.density (pg / um³) = pg
    // cell_density 1.04 pg/um³
    // cell volume ~ 2500 (um³) * ~ 1.04 pg / (um³) 
    // gDW cell (cell.volume.total * cell_density) = mass.total 
    // 2500 (um³) * 1.04 (pg/um³) =~ 2600 (pg) ~= 2.3 ng
    // cell.mass.solid = cell.mass.total * (1-fluid_frac) = 780 pg
 
    
    dFBAIntracellular* dfba_model = static_cast<dFBAIntracellular*>(phenotype.intracellular);
    double current_volume = phenotype.volume.total;
    double Vmax_scale = current_volume / dfba_model->reference_volume;

    std::vector<double> density_vector = pCell->nearest_density_vector(); 

    // unitless ~30%
    float solid_fraction = 1 - phenotype.volume.fluid_fraction;     
    // um³ = um³
    float solid_volume = current_volume * solid_fraction;   
    // um³ * pg / um³ = pg
    float cell_dry_weight = solid_volume * dfba_model->cell_density ;     
    // r = (3V / 4π))^1/3 (um)
    float radius = cbrt( (3./4. * PI * current_volume) ); 
    // um * um = um²
    float cell_surface = 4 * PI * pow(radius, 2);

    // Setp 1 - Update exchange fluxes lower
    //  Standard FBA fluxes units:  mmol /  gDW cell / h
    //  Km: mM = mmol/L
    //  Vmax: mM/h
    //  sdt: surf_density_of_transport
    //  kcat: catalytic constant
    //  Vmax = kcat * sdt * cell_surface

    map<std::string, ExchangeFluxData>::iterator it;
    for(it = dfba_model->substrate_exchanges.begin(); it != dfba_model->substrate_exchanges.end(); it++)
    {
        std::string substrate_name = it->first;
        ExchangeFluxData ex_strut = it->second;

        // geting the amount of substrate
        double substrate_conc = density_vector[ex_strut.density_index];;
        // scaling Vmax ased on cell volume
        double Vmax = ex_strut.Vmax.value * Vmax_scale;
        double Km = ex_strut.Km.value;
        
        // using irreversible Michaelis Menten kinetics to estimate the flux bound
        double max_rate = (Vmax * substrate_conc) / (Km + substrate_conc); // should be calculated from density
        // Change sign to use as lower bound of the exchange flux
        double exchange_flux_lb = -1 * max_rate;
        // Updateing the lower bound of the corresponding exchange flux
        dfba_model->model.setReactionLowerBound(ex_strut.fba_flux_id, exchange_flux_lb);    
    }
}

void update_dfba_outputs(PhysiCell::Cell* pCell, PhysiCell::Phenotype& phenotype, double dt )
{
    // STEP 3 - Update cell volumen using growth rate (first rscale growth rate to 1/min)
    dFBAIntracellular* dfba_model = static_cast<dFBAIntracellular*>(phenotype.intracellular);
    float growth_rate = dfba_model->model.getObjectiveValue();
    
    growth_rate = growth_rate * hours_to_minutes; // growth_rate 1/h -> 1/min
    
    // V(t+1) = V(t) + V(t) * mu * dt = V(t) * (1 + u * dt) 
    float volume_increase_ratio = 1 + growth_rate * dt;
    phenotype.volume.multiply_by_ratio( volume_increase_ratio );
    phenotype.geometry.update(pCell, phenotype, dt);

    float solid_fraction = 1 - phenotype.volume.fluid_fraction;     
    // um³ = um³
    float solid_volume = phenotype.volume.total * solid_fraction;   
    float cell_dry_weight = solid_volume * dfba_model->cell_density ;  

    
    // STEPS 4-5 - Update net_export_rates for the different densities
    map<std::string, ExchangeFluxData>::iterator it;
    for(it = dfba_model->substrate_exchanges.begin(); it != dfba_model->substrate_exchanges.end(); it++)
    {
        // Retrive the exchange flux and its corresponding density
        std::string substrate_name = it->first;
        ExchangeFluxData ex_strut = it->second;
        
        int density_index = ex_strut.density_index;
        std::string fba_flux_id = ex_strut.fba_flux_id;
        
        dFBAReaction* exchange_flux = dfba_model->model.getReaction(fba_flux_id);
        double flux_value = exchange_flux->getFluxValue(); // mmol/gDW/h

        // Rescaling FBA exchanges flux into net_export_rates
        // Net export rates are expressed in substance/time
        // mmol/gDW/h --> mmol/min
        // flux_value / 60 * cell_dry_weight  = pico*mmol/min = fmol/min
        // fmol/min * 10e3 =  pmol/min
        double net_export_rate = flux_value * hours_to_minutes * cell_dry_weight * 1e3;

        phenotype.secretion.net_export_rates[density_index] = net_export_rate;
        
        if (default_microenvironment_options.track_internalized_substrates_in_each_agent)
        {
            phenotype.molecular.internalized_total_substrates[density_index] = 0;
        }
    }

    return;
}

void dFBAIntracellular::save_dFBA(std::string path, std::string index) 
{
	
}

}