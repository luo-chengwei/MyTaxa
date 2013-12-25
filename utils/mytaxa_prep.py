#! /usr/bin/env python

__author__ = 'Chengwei Luo (luo.chengwei@gatech.edu)'
__version__ = '0.0.1'
__date__ = 'November 2013'


"""
This script is part of the MyTaxa package and it's under GNU License v3.0

Copyright(c) 2013 Chengwei Luo (luo.chengwei@gatech.edu), Konstantinidis Laboratory,
			Georgia Institute of Technology.

	This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>

https://github.com/luo-chengwei/MyTaxa

for help, type:
python mytaxa_prep.py --help
"""

USAGE = \
"""Usage: %prog <required_parameters> [options]

Add --help to see a full list of required and optional
arguments to run this script.

If you use MyTaxa in your work, please cite it as:
<MyTaxa citation here>

Copyright: Chengwei Luo, Konstantinidis Lab, Georgia Institute of Technology, 2013
"""

import sys
import os
import re
import glob
import re
from optparse import OptionParser, OptionGroup
from subprocess import PIPE, Popen, call
from Bio import SeqIO

def main(argv = sys.argv[1:]):

	parser = OptionParser(usage = USAGE, version="Version: " + __version__)
	
	# Required arguments
	requiredOptions = OptionGroup(parser, "Required options",
								"These options are required to run BinGeR, and may be supplied in any order.")
	
	requiredOptions.add_option("-i", "--infile", type = "string", metavar = "FILE",
							help = "Input file with query sequences in multi-fasta format.")

	requiredOptions.add_option("-p", "--prefix", type = "string", metavar = "STRING",
							help = "Prefix for output files produced by this script.\
									The files are: <prefix>.gff, <prefix>.faa, and <prefix>.ffn")

	parser.add_option_group(requiredOptions)
	
	# Optional arguments that need to be supplied if not the same as default
	
	optOptions = OptionGroup(parser, "Optional parameters",
						"There options are optional, and may be supplied in any order.")
	
	optOptions.add_option("--prodigal", type = "string", default = "prodigal", metavar = "STRING",
							help = "Location of prodigal executable (Hyatt D, et al, Bioinformatics, 2010).")
	
	parser.add_option_group(optOptions)
	
	(options, args) = parser.parse_args(argv)
	
	if options.infile is None:
		parser.error("An infile is required!")
		exit(0)
		
	if options.prefix is None:
		parser.error("An output prefix is required to supply!")
		exit(0)
	
	# test prodigal
	prodigal_test = Popen(options.prodigal, shell=True, stdout=PIPE, stderr=PIPE).stderr.read()
	if not prodigal_test or prodigal_test.count('not found') ==1:
		sys.stderr.write("FATAL: Prodigal not found in path!\n")
		exit(0)
				
	tempfaa = options.prefix + '.faa.tmp'
	tempffn = options.prefix + '.ffn.tmp'
	faa = options.prefix + '.faa'
	ffn = options.prefix + '.ffn'
	gff = options.prefix + '.gff'
	
	# run prodigal first
#	call([options.prodigal, "-a", tempfaa, "-d", tempffn, "-f", "gff",\
#		 "-o", gff, "-p", "meta", "-i", options.infile],
#		stdout = PIPE, stderr = PIPE)
	
	ofh = open(faa, 'w')	
	for record in SeqIO.parse(tempfaa, 'fasta'):
		tag, start, end, strand = re.search('(.+)\_\d+\s\#\s(\d+)\s\#\s(\d+)\s\#\s(\-?\d+)\s\#\s', record.description).group(1, 2, 3, 4)
		new_name = tag + '|' + start + '-' + end + '|' + strand
		ofh.write('>%s\n%s\n' % (new_name, record.seq))
	ofh.close()
	
	ofh = open(ffn, 'w')	
	for record in SeqIO.parse(tempffn, 'fasta'):
		tag, start, end, strand = re.search('(.+)\_\d+\s\#\s(\d+)\s\#\s(\d+)\s\#\s(\-?\d+)\s\#\s', record.description).group(1, 2, 3, 4)
		new_name = tag + '|' + start + '-' + end + '|' + strand
		ofh.write('>%s\n%s\n' % (new_name, record.seq))
	ofh.close()
	
	os.remove(tempfaa)
	os.remove(tempffn)
	
	
if __name__ == '__main__':
	main()
	
		
	