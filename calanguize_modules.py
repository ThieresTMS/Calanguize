import re
import os
from Bio import SeqIO
from Bio.Seq import Seq


def read_config_file(fh_project, fh_code):
    parameters = {}
    parameters['infile'] = read_line("infile", fh_project)
    parameters['outdir'] = read_line("outdir", fh_project)
    parameters['busco_database'] = read_line("busco_database", fh_project)
    parameters['assembly_summary'] = read_line("assembly_summary", fh_project)
    parameters['extract_mode'] = read_line("extract_mode", fh_project)
    parameters['python'] = read_line("python", fh_code)
    parameters['busco'] = read_line("busco", fh_code)

    return parameters


def read_line(parameter, file):
    for line in file:
        line = line.rstrip()
        match = re.search("^%s" % parameter, line)
        if match:
            line = re.sub('^\s*%s\s*=\s*' % parameter, '', line)
            line = re.sub('#.*$', '', line)
            line = re.sub('\s*$', '', line)
            return line


def check_parameters(parameters, code_config):
    if parameters['infile'] == "":
        exit("No path to the species files was specified in your project's configuration file, please fill the parameter 'infile'.")
    if not os.path.isfile(parameters['infile']):
        exit("The path to your project's species files isn't a valid file, please check if the path in 'infile' is correct: %s" % parameters['infile'])
    if not os.access(parameters['infile'], os.R_OK):
        exit("You don't have permission to read in your project's species file: %s, please redefine your permissions." % parameters['infile'])

    if parameters['outdir'] == "":
        exit("No path to outdir was specified in project_config, please open this file and fill the parameter 'outdir'")

    if parameters['assembly_summary'] == "":
        exit("No path to assembly summary file was specified in your project's configuration file, please fill the parameter 'assembly_summary'.")
    if not os.path.isfile(parameters['assembly_summary']):
        exit("The path to your assembly_summary file isn't a valid file, please check if the path in 'assembly_summary' is correct: %s" % parameters['assembly_summary'])
    if not os.access(parameters['infile'], os.R_OK):
        exit("You don't have permission to read assembly summary file, please redefine your permissions")

    if parameters['python'] == "":
        exit("No path to python was specified in code_config at %s, please open this file and fill the parameter 'python'" % code_config)
    if not os.path.isfile(parameters['python']):
        exit("The executable of python wasn't found in the specified path, please check if the path is correct: %s" % parameters['python'])
    if not os.access(parameters['python'], os.R_OK):
        exit("You don't have permission to execute the python file specified at code_config, please check permissions or replace the file")

    if parameters['busco'] == "":
        exit("No path to busco was specified in code_config at %s, please open this file and fill the parameter 'busco'" % code_config)
    if not os.path.isfile(parameters['busco']):
        exit("The executable of busco wasn't found in the specified path, please check if the path is correct: %s" % parameters['busco'])
    if not os.access(parameters['busco'], os.R_OK):
        exit("You don't have permission to execute the busco file specified at code_config, please check permissions or replace the file")

    if parameters['busco_database'] == "":
        exit("No path to busco database was specified in code_config at %s, please open this file and fill the parameter 'busco_dabase'")
    if not os.path.isdir(parameters['busco_database']):
        exit("The path to busco database wasn't found in the specified path, please check if the path is correct: %s" % parameters['busco_database'])

    if parameters['extract_mode'] == "":
        exit("Extract mode not setted. Please fill the parameter 'extract_mode'")
    if (parameters['extract_mode'] != "protein") and (parameters['extract_mode'] != "orfs"):
        exit("Invalid parameter %s at 'extract_mode'. Please set 'protein' or 'orfs'" % parameters['extract_mode'])


def get_genomes(parameters, species):
    url = {}
    species_in_summary = []
    genomes_out = parameters['outdir']+"/genomes"
    if not os.path.isdir(genomes_out):
        os.mkdir(genomes_out)
    with open(parameters['assembly_summary'], "r") as assembly_summary:
        for line in assembly_summary:
            if line.startswith("#"):
                continue
            line_filds = line.split("\t")
            url[line_filds[7]] = line_filds[19]
            species_in_summary.append(line_filds[7])
    assembly_summary.close()
    ids_find = []
    ids_not_find = []
    for name in species:
        if name in species_in_summary:
            ids_find.append(name)
        else:
            ids_not_find.append(name)
    if ids_not_find:
        ids_err = parameters['outdir']+"/species_not_foud.txt"
        species_not_find = open(ids_err, "a")
        for name in ids_not_find:
            print("I did not find %s in assembly summary, please check if specie name is correct, or if it has a sub-specie name." % name)
            species_not_find.write("%s\n" % name)
        species_not_find.close()
    final_ids = []
    for name in ids_find:
        id = re.sub('\s', '_', name)
        final_ids.append(id)
        path_dowload = genomes_out + "/" + id + ".gbff.gz"
        url_filds = url[name].split("/")
        url[name] = url[name]+"/"+url_filds[-1]+"_rna.gbff.gz"
        print("Downloading genome", name)
        os.system('wget -c -O %s %s' % (path_dowload, url[name]))
        os.system('gzip -d %s' % path_dowload)
    return final_ids


def get_orfs(parameters, id):
    path_to_file_in = parameters['outdir']+"/genomes/"+id+".gbff"
    path_to_file_out = parameters['outdir']+"/orfs/"+id+".fasta"
    out = open(path_to_file_out, "a")
    for seq_record in SeqIO.parse(path_to_file_in, "genbank"):
        seq = ""
        gene_name = ""
        for feature in seq_record.features:
            if 'pseudo' in feature.qualifiers:
                continue
            if feature.type == "CDS":
                seq = feature.location.extract(seq_record).seq
                gene_name = feature.qualifiers['gene'][0]
            else:
                continue
        flag = check_cds(seq)
        if flag is False:
            continue
        out.write(">gene:%s|protein_id:%s\n%s\n" % (gene_name, seq_record.id, seq))
    out.close()


def get_proteins(parameters, id):
    path_to_file_in = parameters['outdir'] + "/genomes/" + id + ".gbff"
    path_to_file_out = parameters['outdir'] + "/proteins/" + id + "_proteins.fasta"
    out = open(path_to_file_out, "a")
    for seq_record in SeqIO.parse(path_to_file_in, "genbank"):
        seq = ""
        gene_name = ""
        for feature in seq_record.features:
            if 'pseudo' in feature.qualifiers:
                continue
            if feature.type == "CDS":
                seq = feature.qualifiers['translation'][0]
                gene_name = feature.qualifiers['gene'][0]
            else:
                continue
        if seq != "":
            out.write(">gene:%s|protein_id:%s\n%s\n" % (gene_name, seq_record.id, seq))
    out.close()


def check_cds(seq):
    start_codon = seq[0:3]
    stop_codon = seq[-3:]
    valid_start_codon = ['ATG', 'GTG']
    if start_codon not in valid_start_codon:
        return False
    valid_stop_codon = ['TAA', 'TAG', 'TGA']
    if stop_codon not in valid_stop_codon:
        return False
    if len(seq) % 3 != 0:
        return False
    if check_code(seq) is False:
        return False
    if len(seq) < 100:
        return False
    return True


def check_code(seq):
    code = "ACTG"
    for base in seq:
        if base not in code:
            return False
    return True


def get_longest(parameters, id):
    gene_data = {}
    seq_in = ""
    out = ""
    if parameters['extract_mode'] == "orfs":
        seq_in = parameters['outdir'] + "/orfs/" + id + ".fasta"
        seq_out = parameters['outdir'] + "/longest_orfs/" + id + "_longest_orfs.fasta"
        out = open(seq_out, "a")
    if parameters['extract_mode'] == "protein":
        seq_in = parameters['outdir'] + "/proteins" + id + "_proteins.fasta"
        seq_out = parameters['outdir'] + "/longest_proteins/" + id + "_longest_proteins.fasta"
        out = open(seq_out, "a")
    for seq_record in SeqIO.parse(seq_in, "fasta"):
        ids_filds = seq_record.id.split("|")
        gene_name = re.sub("gene:", "", ids_filds[0])
        protein_id = re.sub("protein_id:", "", ids_filds[1])
        if gene_name in gene_data:
            actual_length = len(gene_data[gene_name]["sequence"])
            new_length = len(seq_record.seq)
            if new_length > actual_length:
                gene_data[gene_name]["sequence"] = seq_record.seq
                gene_data[gene_name]["protein_id"] = protein_id
        else:
            gene_data[gene_name] = {}
            gene_data[gene_name]["sequence"] = seq_record.seq
            gene_data[gene_name]["protein_id"] = protein_id
    for key in gene_data:
        out.write(">%s\n%s\n" % (gene_data[key]["protein_id"], gene_data[key]["sequence"]))
    out.close()


def translation(parameters, id):
    file_in = parameters['outdir'] + "/longest_orfs/" + id + "_longest_orfs.fasta"
    file_out = parameters['outdir'] + "/longest_proteins/" + id + ".longest_proteins.fasta"
    out = open(file_out, 'a')
    for seq_records in SeqIO.parse(file_in, "fasta"):
        seq = Seq(seq_records.seq)
        protein = seq.translate()
        out.write(">%s\n%s\n" % (seq_records.id, protein))
    out.close()

def run_busco(parameters, id):
    comp_dup = []
    input = parameters['outdir'] + "/longest_proteins/" + id + ".longest_proteins.fasta"
    busco_outdir = parameters['outdir'] + "/busco"
    os.system('%s %s -i %s -o %s -l %s -m prot -c 50  -f' % (parameters['python'], parameters['busco'], input, id, parameters['busco_database']))
    os.system('mv -u run_%s %s' % (id, busco_outdir))
    busco_file = parameters['outdir'] + "/busco/run_" + id + "/short_summary_" + id + ".txt"
    with open(busco_file, "r") as busco:
        for line in busco:
            if re.search("C:", line):
                comp = line.split('%')
                completeness_data  = comp[0].split(':')
                duplicate_data = comp[2].split(':')
                print("Genome %s completeness is %s and has %s duplicate" % (id, completeness_data[1], duplicate_data[1]))

