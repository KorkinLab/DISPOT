#!/usr/local/bin/python
import pickle 
import pandas as pd
import numpy as np
import argparse
import sys
from copy import deepcopy

def normalize_data(domain_dict):
	temperature = 0
	normalization_constant = 0
	original_domain_dict = deepcopy(domain_dict)

	# Dimeric potentials
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


def load_sf_data():
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
		return np.nan, np.nan

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
	parser.add_argument('-s', '--sf', action='store_true', required=False,
							help='Changes default domains from \'Family\' (fa) level of SCOP hierarchy to \'Superfamily\' (sf) level')
	#parser.add_argument('-l', '--list', action='store_true', required=False, 
	#						help='Indicates that a list of domains would be passed. \
	#						Computes statistical potentials for all possible pairs of listed domains.')
	parser.add_argument('-m', '--max', action='store_true', required=False,
							help='For a selected domain produces the highest value of statistical potential and a SCOP \
                             ID for corresponding interaction domain partner.')
	parser.add_argument('-d', '--domains', type=str, nargs='+', required=False,
							help='List of the domains to be considered. \
							Domains should be separated by space. \
							Computes statistical potentials for all possible pairs of listed domains. \
                            If only one domain is provided then statistical potential between the same domains will be returned.')
	parser.add_argument('-i', '--input_file', type=str, nargs=1, required=False,
							help='Input tab-separated file with lines in format <domain1>	<domain2> for \
							pairwise statistic. Or in a single line format <domain1>	<domain2>	<domain3> 	... \
							for maximal potential statistic. It is an alternative to the --domains input option.')
	parser.add_argument('-o', '--output', type=str, nargs=1, required=False,
							help='Path to the output file. By default the entire output is printed on the console (stdout). \
							Attention - all directories in the file path should exist.')
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


def run(args):
	data = None
	if args.sf:
		data = load_sf_data()
	else:
		data = load_data()

	data, monomeric_data = normalize_data(data)
	monomeric_domains = None
	domain_pairs = None
	max_domains = None
	out = None

	if args.output is not None:
		out = open(args.output[0], 'w')

	if args.domains is not None:
		if args.mono:
			monomeric_domains = args.domains
		elif args.max:
			max_domains = args.domains
		else:
			domain_pairs = []
			if len(args.domains) < 2:
				raise Exception('For pairwise statistical potentials at least 2 arguments should be passed')
			for i in range(len(args.domains)):
				for j in np.arange(i, len(args.domains)):
					pair = (args.domains[i], args.domains[j])
					reverse_pair = (pair[1], pair[0])
					if pair not in domain_pairs and reverse_pair not in domain_pairs:
						domain_pairs.append(pair)

	if args.input_file is not None:
		domains = pd.read_csv(args.input_file[0], sep='\t')
		domains = np.squeeze(domains.as_matrix())
		domains = domains.astype(str)
		if args.mono:
			if monomeric_domains is None:
				monomeric_domains.append(domains)
		if args.max:
			if max_domains is None:
				max_domains = []
			max_domains.append(domains)
		else:
			if np.shape(domains)[1] < 2:
				raise Exception('For pairwise statistical potentials at least 2 arguments should be passed')
			if domain_pairs is None:
				domain_pairs = []
			for i in range(np.shape(domains)[0]):
				pair = (domains[i][0], domains[i][1])
				reverse_pair = (pair[1], pair[0])
				if pair not in domain_pairs and reverse_pair not in domain_pairs:
					domain_pairs.append(pair)

	if monomeric_domains is not None:
		process_monomeric_potentials(monomeric_data, monomeric_domains)
	if max_domains is not None:
		process_max_potentials(data, max_domains)
	if domain_pairs is not None:
		process_domain_pairs(data, domain_pairs)
	if max_domains is None and domain_pairs is None and monomeric_domains is None:
	 		raise Exception('No data passed. Please use --domain or --input_file option')

if __name__ == "__main__":
	args = parse_arguments(sys.argv[1:])
	run(args)