import argparse
import sys
import os

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
	args = parser.parse_args(arguments)
	parser.description = "This is a docker version of DISPOT software. For correct work it requires installed Docker and configured docker group. \
							Information about installing and configuring Docker: https://docs.docker.com/install/. Select your operation system (OS) \
							and proceed with instructions."
	if len(arguments) == 0:
		parser.print_help(sys.stderr)
		sys.exit(1)
	return args


def remove_output(input):
	for i in range(len(input)):
		if input[i] == "-o" or input[i] == "--output":
			if i == (len(input) - 2):
				return input[:i], input[i+1]
			else:
				return input[:i] + input[i+2:], input[i+1]
	return input, None

if __name__ == "__main__":
	docker_image = 'narykov/dispot:basic'
	args = parse_arguments(sys.argv[1:])
	new_args, output = remove_output(sys.argv[1:])
	print(sys.argv[1:])
	print(new_args)
	if args.output is not None:
		os.system('docker run {} {} > {}'.format(docker_image, ' '.join(new_args), output))
	else:
		os.system('docker run {} {}'.format(docker_image, ' '.join(new_args)))