                       README for Calanguize
             (last updated 11/30/2020 - mm/dd/yyyy format)

AUTHORS
-=-=-=-

Thieres Tayroni Martins da Silva (thierestayroni@gmail.com)


1 - DESCRIPTION
-=-=-=-=-=-=-=-

Calanguize is software that prepares your data to be used in Calango.


2 - HOW TO USE - IN BRIEF
-=-=-=-=-=-=-=-=-=-=-=-=-

2.1 - INSTALL 
-=-=-=-=-=-=-

Download calanguize's dependencies detailed in the INSTALL file and set the
parameters in 'komodize_config' file .


2.2 - PREPARING INPUT FILES 
-=-=-=-=-=-=-=-=-=-=-=-=-=-

The input file for Komodize is very simple, consists of a list containing the scientific names of the species to be analyzed.
These names must follow the binomial nomenclature rule, where the first name must begin with capital letter and the second with a lowercase letter separated by space.
If the species has a subspecies, the full name of the subspecies should be written (e.g. Canis lupus familiaris or Gorilla gorilla gorilla)

 - Configure a project configuration file.

Copy the project_config file and enter the name of your project "cp project_config $project_name" and set the parameters.

infile =     #path to file that contains species names to be analyzed

outdir =     #path to output directory

busco_database = #path to busco datasets

max_processors =    #maximum number of cpus to be used

extract_mode =  #sequences to be extracted from the GBFF file, protein or orfs;
    

2.3 - RUNNING Calanguize 
-=-=-=-=-=-=-=-=-=-=-    

Execute Calanuize, using as input the project configuration file you generated
in the previous step:

  python3.8 calanguize_genomes.py code_config project_config
