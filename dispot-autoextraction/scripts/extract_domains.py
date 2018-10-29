import os
import sys
import pickle
import subprocess
import numpy as np
import pandas as pd
from Bio import SeqIO


def make_fasta(seq, fname=None):
	if fname is None:
		fname = 'seq_fasta'
		fpath = 'tmp/seq_fasta'
	else:
		fpath = 'tmp/{}'.format(fname)
	f = open(fpath, 'w')
	f.write('>SEQ\n')
	f.write(seq)
	f.close()
	return fpath, fname

def parse_fasta(fasta_files, ignore_decimals=True):
	fasta_records = {}
	for f in fasta_files:
		record_dict = SeqIO.to_dict(SeqIO.parse(f, "fasta"))
		if ignore_decimals:
			for key in record_dict.keys():
				fasta_records[key.split('.')[0]] = record_dict[key]
		else:
			fasta_records = merge_two_dicts(fasta_records, record_dict)
	return fasta_records

def merge_two_dicts(x, y):
    z = x.copy()   # start with x's keys and values
    z.update(y)    # modifies z with y's keys and values & returns None
    return z

def extract_domains_superfamily(seq, name):
	fasta, fasta_name = make_fasta(seq)
	cur_dir = os.getcwd()
	os.system('cp {} ./superfamily/.'.format(fasta))
	os.chdir('./superfamily')
	os.system('perl superfamily.pl {}'.format(fasta_name))
	os.system('rm ./scratch/*')
	os.chdir(cur_dir)
	os.system('mv ./superfamily/{}.ass tmp/{}.ass'.format(fasta_name, fasta_name))
	os.system('python ./scripts/parse_domain_information.py ./tmp/{}.ass'.format(fasta_name))
	domains = pickle.load(open('tmp/parsed.pickle', 'r'))
	binding_sites = pickle.load(open('tmp/binding_sites.pickle', 'r'))
	return domains, binding_sites


def load_stat_pot_data():
	data = pd.read_csv('./data/curated_statistical_potentials.csv', header=None)
	output = {}
	for i in data.index:
		pair, potential = data.loc[i, data.columns]
		domain1, domain2 = pair.split('-')
		if not output.has_key(domain1):
			output[domain1] = {}
		output[domain1][domain2] = potential
		if not output.has_key(domain2):
			output[domain2] = {}
			output[domain2][domain1] = potential
	return output



if __name__ == "__main__":
	path = None
	if len(sys.argv) != 2:
		print("Error. The program takes in the following input: <path to a directory with FASTA files>")
		raise Exception("Wrong input")
	else: 
		path = sys.argv[1]
	fasta_files = os.listdir(path)
	fasta_files = [os.path.join(path, x) for x in fasta_files if os.path.isfile(os.path.join(path, x))]
	print(fasta_files)
	fasta = parse_fasta(fasta_files)
	for key in fasta.keys():
		seq = str(fasta[key].seq)
		name = key.split('|')[0]
		name = name.split(':')[0] + '_' + name.split(':')[1]
		domains , binding_sites = extract_domains_superfamily(seq, name)
		print("Name: {}".format(name))
		print("BInding Sites: {}".format(binding_sites))
		pickle.dump(binding_sites, open('tmp/{}.bs'.format(name), 'w+'))