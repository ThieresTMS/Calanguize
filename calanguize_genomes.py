import sys
import calanguize_modules
import os

project_config = sys.argv[2]
code_config = sys.argv[1]

try:
    fh_project_config = open(project_config, "r")
except OSError:
    print("Couldn't open the project configuration file, you may have to chek if the name or permission of the file are correct")
    sys.exit()
try:
    fh_code_config = open(code_config, "r")
except OSError:
    print("Couldn't open the config file, you may have to chek if the name or permission of the file are correct")
    sys.exit()

parameters = calanguize_modules.read_config_file(fh_project_config, fh_code_config)
fh_project_config.close()
fh_code_config.close()

calanguize_modules.check_parameters(parameters, code_config)

species = []
with open(parameters['infile'], "r") as specie_list:
    for line in specie_list:
        line = line.rstrip()
        species.append(line)
specie_list.close()

print("Welcome to Calanguize")
if os.path.isdir(parameters['outdir']):
    exit("Outdir '%s' already exists, please change the parameter 'outdir'," % parameters['outdir'])
else:
    os.mkdir(parameters['outdir'])

ids = calanguize_modules.get_genomes(parameters, species)


if parameters['extract_mode'] == "orfs":
    orfs_out = parameters['outdir'] + "/orfs"
    if not os.path.isdir(orfs_out):
        os.mkdir(orfs_out)
    for id in ids:
        print("Extracting ORFS from %s" % id)
        calanguize_modules.get_orfs(parameters, id)
    orfs_longest_out = parameters['outdir'] + "/longest_orfs"
    if not os.path.isdir(orfs_longest_out):
        os.mkdir(orfs_longest_out)
    for id in ids:
        print("Summarizing ORFS from %s per locus" % id)
        calanguize_modules.get_longest(parameters, id)
    translate = parameters['outdir'] + "/longest_proteins"
    if not os.path.isdir(translate):
        os.mkdir(translate)
    for id in ids:
        print("Translating ORFs from %s" % id)
        calanguize_modules.translation(parameters, id)

if parameters['extract_mode'] == "protein":
    proteins_out = parameters['outdir'] + "/proteins"
    if not os.path.isdir(proteins_out):
        os.mkdir(proteins_out)
    for id in ids:
        print("Extracting proteins from %s" % id)
        calanguize_modules.get_proteins(parameters, id)
    proteins_longest_out = parameters['outdir'] + "/longest_proteins"
    if not os.path.isdir(proteins_longest_out):
        os.mkdir(proteins_longest_out)
    for id in ids:
        print("Summarizing proteins from %s per locus" % id)
        calanguize_modules.get_longest(parameters, id)


busco_out = parameters['outdir'] + "/busco"
if not os.path.isdir(busco_out):
    os.mkdir(busco_out)
for id in ids:
    print("Running busco for %s" % id)
    calanguize_modules.run_busco(parameters, id)

print("Finish! Closing Komodize, see you around")
