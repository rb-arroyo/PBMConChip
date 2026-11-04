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
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
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

#include "./custom.h"

// Global variable
std::vector<std::vector<double>> t_cell_positions;

void create_cell_types( void )
{
	// set the random seed 
	if (parameters.ints.find_index("random_seed") != -1)
	{
		SeedRandom(parameters.ints("random_seed"));
	}
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 
	
	initialize_default_cell_definition(); 
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	cell_defaults.functions.contact_function = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 

	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	build_cell_definitions_maps(); 

	/*
	   This intializes cell signal and response dictionaries 
	*/

	setup_signal_behavior_dictionaries(); 

/*
       Cell rule definitions 
	*/

	setup_cell_rules(); 
	

	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 
	
	cell_defaults.functions.update_phenotype = phenotype_function; 
	cell_defaults.functions.custom_cell_rule = custom_function; 
	cell_defaults.functions.contact_function = contact_function; 
	
	Cell_Definition* pCD = find_cell_definition( "PBMC" );
	double mr = pCD->custom_data["multiplication_ratio_value"];
        pCD-> phenotype.volume.multiply_by_ratio(mr);
	pCD->functions.update_migration_bias = custom_migration_speed;
	
	
	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	display_cell_definitions( std::cout ); 
	
	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters 
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	
	
	return; 
}

void setup_tissue( void )
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
	
	/*Cell* pCell;
	
	Cell_Definition* pCD = find_cell_definition( "PBMC" );
	pCell = create_cell( *pCD );
	pCell->assign_position( -1.0 , -100.0 , 0.0 );
	std::cout << std::endl;*/
	
	// upload a file with {t,x,y,z} for each cell 
	t_cell_positions = load_file_vector("./config/cells_opt_c4_d1.csv");

	std::cout << std::endl; 
	
	// load cells from your CSV file (if enabled)
	load_cells_from_pugixml(); 	
	
	return; 
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{ return paint_by_number_cell_coloring(pCell); }

void custom_migration_speed( Cell* pCell, Phenotype& phenotype, double dt )
{
	// find my cell definition
	Cell_Definition* pCD = find_cell_definition( pCell->type );
	
	// set migration bias direction to grad(pO2)
	phenotype.motility.migration_bias_direction = {-1,0,0};
	
	// trick: if dead, overwrite with NULL function pointer.
	if( phenotype.death.dead == true )
	{ pCell->functions.update_migration_bias= NULL; }
return;
}

void create_exp_cell(double t_create)
{
	Cell* pCell;
	Cell_Definition* pCD = find_cell_definition( "PBMC" );
	
	for( int n = 0 ; n < t_cell_positions.size()  ; n++ )
	{	
		if (fabs(t_cell_positions[n][0] - t_create) < 0.02)
		{
			std::vector<double> position = {0,0,0}; 
			position[0] = t_cell_positions[n][1];
			position[1] = t_cell_positions[n][2];
			position[2] = t_cell_positions[n][3];

			pCell = create_cell( *pCD ); 
			pCell->assign_position( position );
			
			//Aseguramos que la célula sea móvil
            		pCell->phenotype.motility.is_motile = true; 
            
            		// Establecemos la direccion del vector dirigido (d_bias)
            		pCell->phenotype.motility.migration_bias_direction = {-1,0,0}; 
            
           		 //Forzamos el cálculo inmediato de la velocidad.
           		 //Esto evita que espere al siguiente ciclo de persistencia.(u<dt/Tpers)
            		custom_update_motility_vector(pCell,pCell->phenotype, mechanics_dt );
            		
			//std::cout<< "A cell created " << std::endl;
		}
	}
	
}

void custom_update_motility_vector(Cell* pCell, Phenotype& phenotype , double dt_)
{	if( phenotype.motility.is_motile == false )
	{
		phenotype.motility.motility_vector.assign( 3, 0.0 ); 
		return; 
	}
	std::vector<double> randvec(3,0.0);
		if( phenotype.motility.restrict_to_2D == true )
		{ randvec = UniformOnUnitCircle(); }
		else
		{ randvec = UniformOnUnitSphere(); }

		// if the update_bias_vector function is set, use it  
		if( pCell->functions.update_migration_bias )
		{
			pCell->functions.update_migration_bias( pCell,phenotype,dt_ ); 
		}
		
		phenotype.motility.motility_vector = phenotype.motility.migration_bias_direction; // motiltiy = bias_vector
		phenotype.motility.motility_vector *= phenotype.motility.migration_bias; // motility = bias*bias_vector 
		
		double one_minus_bias = 1.0 - phenotype.motility.migration_bias; 
		
		axpy( &(phenotype.motility.motility_vector), one_minus_bias, randvec ); // motility = (1-bias)*randvec + bias*bias_vector
		
		normalize( &(phenotype.motility.motility_vector) ); 
		
		phenotype.motility.motility_vector *= phenotype.motility.migration_speed; 

return;}

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{ return; }

void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{ return; } 

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{ return; }

std::vector<std::vector<double>> load_file_vector(std::string filename)
{
    std::vector<std::vector<double>> t_cell_positions;
    std::ifstream file(filename);

    std::string line;
    
     if (!file.is_open())
    {
        std::cerr << "ERROR: cannot open file " << filename << std::endl;
        return {};
    }
   
    while (std::getline(file, line))
    {
        std::stringstream ss(line);
        std::string t_str, x_str, y_str, z_str;

        if (std::getline(ss, t_str, ',') && std::getline(ss, x_str, ',') && std::getline(ss, y_str, ',') && std::getline(ss, z_str, ','))
        {
            double t = std::stod(t_str);
            double x = std::stod(x_str);
            double y = std::stod(y_str);
            double z = std::stod(z_str);
            t_cell_positions.push_back({t, x, y, z});
        }
    }

    file.close();
    
    return t_cell_positions;
} 
