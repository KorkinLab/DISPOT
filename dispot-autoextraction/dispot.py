import pickle 
import pandas as pd
import numpy as np
import argparse
import sys
import os
from copy import deepcopy

def normalize_data(domain_dict):
	temperature = 0
	normalization_constant = 0
	original_domain_dict = deepcopy(domain_dict)

	# DImeric potentials
	values = []
	for key1 in domain_dict.keys():
		for key2 in domain_dict[key1].keys():
			value = float(domain_dict[key1][key2])
			values.append(value)

	temperature = np.mean(values)
	for value in values:
		normalization_constant += np.log(float(value)/temperature)

	for key1 in domain_dict.keys():
		for key2 in domain_dict[key1].keys():
			domain_dict[key1][key2] = np.log(domain_dict[key1][key2]/temperature)

	fm = open('monodimeric_stat_potentials.tsv', 'w')
	fh = open('heterodimeric_stat_poteintials.tsv', 'w')	

	triplets = []
	for key1 in domain_dict.keys():
		for key2 in domain_dict[key1].keys():
			triplet = (key1, key2, domain_dict[key1][key2])
			reversed_triplet = (triplet[1], triplet[0], triplet[2])
			if triplet not in triplets and reversed_triplet not in triplets:
				triplets.append(triplet)
	triplets = sorted(triplets, key = lambda t: t[2], reverse = True)

	for triplet in triplets:
		if triplet[0] == triplet[1]:
			fm.write('{}\t{}\t{}\n'.format(triplet[0], triplet[1], triplet[2]))
		else:
			fh.write('{}\t{}\t{}\n'.format(triplet[0], triplet[1], triplet[2]))

	# Monomeric potentials
	monomeric_domain_dict = {}
	values = []
	for key1 in original_domain_dict.keys():
		N = 0
		for key2 in original_domain_dict[key1].keys():
			N += float(original_domain_dict[key1][key2])
		monomeric_domain_dict[key1] = N 
		values.append(N)

	temperature = np.mean(values)
	for value in values:
		normalization_constant += np.log(float(value)/temperature)
	
	for key in monomeric_domain_dict.keys():
		monomeric_domain_dict[key] = np.log(monomeric_domain_dict[key]/temperature)

	return domain_dict, monomeric_domain_dict


def write_pair_with_no_data(data, stat_pot):
	output_path = 'no_interaction_pairs.tsv'
	out = open(output_path, 'w')
	domains = list(data[data.columns[0]]) + list(data[data.columns[1]])
	for i in range(len(domains)):
		for j in np.arange(i, len(domains)):
			if stat_pot.has_key(domains[i]):
				if stat_pot[domains[i]].has_key(domains[j]):
					continue
			out.write('{}\t{}\tnan'.format(domains[i], domains[j]))
			if domains[i] != domains[j]:
				out.write('{}\t{}\tnan'.format(domains[j], domains[i]))


def family2superfamily_dict():
	scope_cla_path = 'db/dir.cla.scop.txt'
	data = pd.read_csv(scope_cla_path, sep='\t', header=None, usecols=[7,8,9,10,11])
	data = data.as_matrix()
	fa2sf = {'-': '-'}
	N, M = np.shape(data)
	for i in range(N):
		fa_idx = None
		for j in range(M):
			if 'fa=' in data[i,j]:
				fa_idx = j
		fa2sf[str(data[i,fa_idx].split('=')[-1].replace('\n',''))] = str(data[i,fa_idx-1].split('=')[-1].replace('\n', ''))
	return fa2sf


def load_ssf_data():
	data = pd.read_csv('db/curated_statistical_potentials.csv', header=None)
	fa2sf = family2superfamily_dict()
	output = {}
	for i in data.index:
		pair, potential = data.loc[i, data.columns]
		domain1, domain2 = pair.split('-')
		domain1 = fa2sf[str(domain1)]
		domain2 = fa2sf[str(domain2)]
		if not output.has_key(domain1):
			output[domain1] = {}
		if not output[domain1].has_key(domain2):
			output[domain1][domain2] = potential
		else:
			output[domain1][domain2] += potential
		if not output.has_key(domain2):
			output[domain2] = {}
		if not output[domain2].has_key(domain1):
			output[domain2][domain1] = potential
		else:
			output[domain2][domain1] += potential
	return output


def load_data():
	data = pd.read_csv('db/curated_statistical_potentials.csv', header=None)
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


def get_potential(data, domain1, domain2):
	domain1 = str(domain1)
	domain2 = str(domain2)
	if data.has_key(domain1):
		if data[domain1].has_key(domain2):
			return data[domain1][domain2]
	return np.nan


def get_maximal_statistical_potential(data, query_domain):
	if query_domain not in data.keys():
		return np.nan

	potentials = {}
	for key in data[query_domain].keys():
		potentials[key] = data[query_domain][key]

	res = max(potentials.iteritems(), key=lambda (k,v): (v,k))
	return res[1], res[0]

def get_monomeric_statistical_potential(data, query_domain):
	query_domain = str(query_domain)
	if data.has_key(query_domain):
		return data[query_domain]
	return np.nan

def parse_arguments(arguments):
	parser = argparse.ArgumentParser()
	#parser.add_argument('-s', '--ssf', action='store_true', required=False,
	#						help='Changes default domains from family (fa) resolution to superfamily (ssf)')
	parser.add_argument('-m', '--max', action='store_true', required=False,
							help='For a selected domain produces the highest value of statistical potential and a SCOP \
                             ID for corresponding interaction domain partner.')
	parser.add_argument('-f', '--fasta_folder', type=str, nargs=1, required=False,
							help='Path to the folder with FASTA files')
	parser.add_argument('-o', '--output_folder', type=str, nargs=1, required=False,
							help='Path to the output folder. By default all output is printed on the console (stdout). \
							Attention - all directories in a file path should exist.')
	parser.add_argument('-n', '--mono', action='store_true', required=False,
							help='Returns monomeric statistical potentials. Input is the same as with the usage of \'--max\' argument.')
	parser.description = "DISPOT is a software which provides statistical potentials for protein domains."
	args = parser.parse_args(arguments)
	if len(arguments) == 0:
		parser.print_help(sys.stderr)
		sys.exit(1)
	return args


def write(string, f=None):
	if f is None:
		print string,
	else:
		f.write(string)

def process_monomeric_potentials(monomeric_data, domains, out=None):
	for domain in domains:
		monomeric_potential = get_monomeric_statistical_potential(monomeric_data, domain)
		write('{}\t{}\n'.format(domain, monomeric_potential), out)

def process_max_potentials(data, domains, out=None):
	for domain in domains:
		max_pot, interactor = get_maximal_statistical_potential(data, domain)
		write('{}\t{}\t{}\n'.format(domain, interactor, max_pot), out)


def process_domain_pairs(data, domain_pairs, out=None):
	for domain_pair in domain_pairs:
		stat_pot = get_potential(data, domain_pair[0], domain_pair[1])
		write('{}\t{}\t{}\n'.format(domain_pair[0], domain_pair[1], stat_pot), out)

def load_domains_bs(path):
	domains_bs = {}
	files = os.listdir(path)
	for fname in files:
		if '.bs' not in fname:
			continue
		structure_name = fname.split('_')[0].lower()
		chain = fname.split('_')[1].split('.')[0].upper()
		bs = pickle.load(open(os.path.join(path, fname), 'r'))
		if len(bs) == 0:
			continue
		if structure_name in domains_bs.keys():
			domains_bs[structure_name][chain] = bs['SEQ']#restructure_bs(bs['SEQ'])
		else:
			domains_bs[structure_name] = {}
	return domains_bs


def generate_domain_pairs(bs_dict):
	domain_pairs_dict = {}
	domains_total = {}
	#print('BS dict: {}'.format(bs_dict))
	for structure in bs_dict.keys():
		domain_pairs_dict[structure] = []
		domains_total[structure] = []
		structure_domains = []
		for chain in bs_dict[structure].keys():
			dom = bs_dict[structure][chain].keys()
			if len(dom) < 1:
				continue
			if dom[0] == '-':
				continue
			structure_domains += dom#bs_dict[structure][chain].keys()
		#print("Structure domains")
		#print(structure_domains)
		for i in range(len(structure_domains)):
			domains_total[structure].append(structure_domains[i])
			for j in np.arange(i+1, len(structure_domains)):
				pair = (structure_domains[i], structure_domains[j])
				reverse_pair = (structure_domains[j], structure_domains[i])
				if pair not in domain_pairs_dict[structure] and reverse_pair not in domain_pairs_dict[structure]:
					#print(pair)
					domain_pairs_dict[structure].append(pair)
	return domain_pairs_dict, domains_total


def run(args):
	bs_output_path = 'tmp'
	data = None
	fa2sf = None
	max_domains = None
	data = load_ssf_data()
	fa2sf = family2superfamily_dict()

	data, monomeric_data = normalize_data(data)
	monomeric_domains = None
	domain_pairs_dict = None
	max_domains = None
	out = None

	if args.output_folder is not None:
		out = args.output_folder[0]
	else:
		out = '.'

	if args.fasta_folder is not None:
		os.system('python scripts/extract_domains.py {}'.format(args.fasta_folder[0]))
		bs_dict = load_domains_bs(bs_output_path)
		domain_pairs_dict, domains_dict = generate_domain_pairs(bs_dict)
		for structure in domain_pairs_dict.keys():
			domain_pairs = domain_pairs_dict[structure]
			domains = domains_dict[structure]
			if args.mono:
				if monomeric_domains is None:
					monomeric_domains = []
				monomeric_domains += domains
			if args.max:
				if max_domains is None:
					max_domains = []
				max_domains += domains
			if domain_pairs is None:
				domain_pairs = []

			print('Cur Dir Path: {}'.format(os.getcwd()))
			if monomeric_domains is not None:
				output_file = os.path.join(out, '{}.mono'.format(structure))
				output = open(output_file, 'w')
				process_monomeric_potentials(monomeric_data, monomeric_domains, output)
				output.close()
			if max_domains is not None:
				output_file = os.path.join(out, '{}.max'.format(structure))
				output = open(output_file, 'w')
				process_max_potentials(data, max_domains, output)
				output.close()
			if domain_pairs is not None:
				output_file = os.path.join(out, "{}.pot".format(structure))
				output = open(output_file, 'w')
				process_domain_pairs(data, domain_pairs, output)
				output.close()
			if max_domains is None and domain_pairs is None and monomeric_data is None:
				raise Exception('No data passed. Please use --input_folder option')


if __name__ == "__main__":
	tmp_dir = 'tmp'
	if len(os.listdir(tmp_dir)) > 0:
		os.system("rm {}/*".format(tmp_dir))
	args = parse_arguments(sys.argv[1:])
	run(args)