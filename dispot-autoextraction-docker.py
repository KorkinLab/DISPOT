import argparse
import sys
import os

def parse_arguments(arguments):
	parser = argparse.ArgumentParser()
	parser.add_argument('-m', '--max', action='store_true', required=False,
							help='For a selected domain produces the highest value of statistical potential and a SCOP \
							 ID for corresponding interaction domain partner.')
	parser.add_argument('-n', '--mono', action='store_true', required=False,
							help='Returns monomeric statistical potentials. Input is the same as with the usage of \'--max\' argument.')
	parser.add_argument('-f', '--fasta_folder', type=str, nargs=1, required=True,
							help='Path to the folder with FASTA files')
	args = parser.parse_args(arguments)
	parser.description = "This is a docker version of DISPOT software. For correct work it requires installed Docker and configured docker group. \
							Information about installing and configuring Docker: https://docs.docker.com/install/. Select your operation system (OS) \
							and proceed with instructions."
	if len(arguments) == 0:
		parser.print_help(sys.stderr)
		sys.exit(1)
	return args


def remove_fasta_folder(input):
	for i in range(len(input)):
		if input[i] == "-f" or input[i] == "--fasta_folder":
			if i == (len(input) - 2):
				return input[:i], input[i+1]
			else:
				return input[:i] + input[i+2:], input[i+1]
	return input, None

if __name__ == "__main__":
	docker_image = 'narykov/dispot:autoextraction'
	args = parse_arguments(sys.argv[1:])
	new_args, fasta_folder = remove_fasta_folder(sys.argv[1:])
	print(sys.argv[1:])
	print(new_args)
	if args.fasta_folder is not None:
		os.system('docker run -v {}:/dispot/data {} {}'.format(fasta_folder, docker_image, ' '.join(new_args)))
	else:
		raise Exception("Folder with fasta files must be provided. Use argument -f or --fasta_folder")