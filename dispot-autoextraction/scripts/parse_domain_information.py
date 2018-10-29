import pickle
import numpy as np
import pandas as pd
import sys
import os


def parse_interpro(files):
	parsed = {}
	binding = {}
	for fname in files:
		f = open(fname)
		for line in f:
			splitted = line.split('\t')
			split_size = len(splitted)
			id = splitted[0]
			domain = splitted[4]
			if id in parsed.keys():
				parsed[id].append(domain)
			else:
				parsed[id] = [domain]
			binding_sites = [splitted[6], splitted[7]]
			if id in binding.keys():
				binding[id][domain] = binding_sites
			else:
				binding[id] = {}
				binding[id][domain] = binding_sites
	pickle.dump(parsed, open('../tmp/parsed.pickle', 'w'))
	pickle.dump(binding, open('../tmp/binding_sites.pickle', 'w'))

def family2superfamily_dict(scope_cla_path):
	print('Cur wd: {}'.format(os.getcwd()))
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

def parse_superfamily(files, fa2sf=None):
	parsed = {}
	binding = {}
	for fname in files:
		f = open(fname)
		for line in f:
			splitted = line.split('\t')
			id = splitted[0]
			print(splitted)
			domain = None
			if fa2sf is not None:
				domain = fa2sf[str(splitted[-1].replace('\n',''))]
			else:
				domain = str(splitted[-1].replace('\n', ''))
			for bs in splitted[2].split(','):
				binding_sites = bs.split('-')
				if id in parsed.keys():
					parsed[id].append(domain)
				else:
					parsed[id] = [domain]
				if id in binding.keys():
					if domain in binding[id].keys():
						binding[id][domain] += binding_sites
					else:
						binding[id][domain] = binding_sites
				else:
					binding[id] = {}
					binding[id][domain] = binding_sites
	pickle.dump(parsed, open('./tmp/parsed.pickle', 'w'))
	pickle.dump(binding, open('./tmp/binding_sites.pickle', 'w'))

if __name__ == "__main__":
	files = []
	if len(sys.argv) > 2:
		files =  sys.argv[1:]
	else:
		files = [sys.argv[1]]
	fa2sf = family2superfamily_dict('./superfamily/dir.cla.scope.2.06-stable.txt')
	parse_superfamily(files, fa2sf)