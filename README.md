DISPOT software provides statistical potentials for domain interactions. 

It is available as standalone Python 2.7 scripts (folders dispot-basic and dispot-autoextraction) and python wrappers over docker container. 

It is possible to use Docker versions without downloading rest of the DISPOT code. The only prerequisite is an installed docker.

### Basic Docker Version

dispot-docker.py is a standalone script that would download image from DockerHub. The only dependency it have is Docker. For correct work it requires installed Docker and configured docker group. Information about installing and configuring Docker: https://docs.docker.com/install/. Select your operation system (OS) and proceed with instructions. Usage is the same as for the Python version:

```
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
                        ```


### Autoextraction Docker Version
