#!/usr/local/bin/python
import os
import sys
import argparse


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('-m', '--max', action='store_true', required=False,
							help='Produces maximal statistical potential for each domain instead of \
							pairwise potentials if listed.')
	parser.add_argument('-n', '--mono', action='store_true', required=False,
							help='Returns monomeric statistical potentials. Input is the same as with the usage of \'--max\' argument.')
	#parser.add_argument('-s', '--ssf', action='store_true', required=False,
	#						help='Changes default domains from family (fa) resolution to superfamily (ssf)')
	args = parser.parse_args(sys.argv[1:])
	flags = '{m}{{n}}'
	if args.max:
		flags = flags.format(m='-m')
	else:
		flags = flags.format(m='')
	if args.mono:
		flags = flags.format(n=' -n')
	else:
		flags = flags.format(n='')
	print(flags)
	os.system('python dispot.py -f ./data -o ./data/results {}'.format(flags))