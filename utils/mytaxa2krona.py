#! /usr/bin/env python

__author__ = 'Chengwei Luo (luo.chengwei@gatech.edu)'
__version__ = '0.0.1'
__date__ = 'November 2013'

USAGE = \
"""
The usage of this script is:
python mytaxa2krona.py <mytaxa_output> <krona_input>

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
"""

import sys
import os
import re
import glob
import re

def main():
	if len(sys.argv) < 3:
		sys.stderr.write(USAGE)
		exit(0)
	
	ifh = open(sys.argv[1], 'r')
	ofh = open(sys.argv[2], 'w')
	
	lib = {}
	
	while 1:
		line = ifh.readline().rstrip('\n')
		if not line:
			break
		tax = ifh.readline().rstrip('\n')
		if tax not in lib:
			lib[tax] = 0
		lib[tax] += 1
	ifh.close()
	
	for tax in lib:
		count = lib[tax]
		ele = tax.split(';')
		col = []
		for i in range(len(ele)):
			c = re.sub('\<.+\>', '', ele[i])
			col.append(c)
		ofh.write(str(count)+'\t'+'\t'.join(col)+'\n')
	
	ofh.close()
			
if __name__ == '__main__':
	main()
