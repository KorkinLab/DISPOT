DISPOT software provides statistical potentials for domain interactions.

DISPOT takes SUPERFAMILY identifiers for both ‘Family’ (fa) or ‘Superfamily’ (sf) level of SCOP hierarchy as an input and produces statistical potential for corresponding pair of domains. Switching between the SCOP levels is implemented in command line option ssf.  One of the possible input options is a command line option domains, which provides a list of space separated SCOP identifiers. Based on this list, the program produces all possible unique pairwise combinations of identifiers and the corresponding statistical potentials. Option max produces the highest value of statistical potential for a selected domain and a SCOP ID for the corresponding interaction domain partner. Option output specifies the output file. If no file path was specified, then program opens a console output prompting a user to input the data. 


### Python version ###
'''
usage: dispot.py [-h] [-s] [-m] [-d DOMAINS [DOMAINS ...]] [-i INPUT_FILE]
                 [-o OUTPUT]

optional arguments:
  -h, --help            show this help message and exit
  -s, --sf              Changes default domains from 'Family' (fa) level of
                        SCOP hierarchy to 'Superfamily' (sf) level
  -m, --max             For a selected domain produces the highest value of
                        statistical potential and a SCOP ID for corresponding
                        interaction domain partner.
  -d DOMAINS [DOMAINS ...], --domains DOMAINS [DOMAINS ...]
                        List of the domains to be considered. Domains should
                        be separated by space. Computes statistical potentials
                        for all possible pairs of listed domains. If only one
                        domain is provided then statistical potential between
                        the same domains will be returned.
  -i INPUT_FILE, --input_file INPUT_FILE
                        Input tab-separated file with lines in format
                        <domain1> <domain2> for pairwise statistic. Or in a
                        single line format <domain1> <domain2> <domain3> ...
                        for maximal potential statistic. It is an alternative
                        to the --domains input option.
  -o OUTPUT, --output OUTPUT
                        Path to the output file. By default the entire output
                        is printed on the console (stdout). Attention - all
                        directories in the file path should exist.
  -n, --mono            Returns monomeric statistical potentials. Input is the
                        same as with the usage of '--max' argument.

'''
