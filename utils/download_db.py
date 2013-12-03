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

"""

USAGE = \
"""Usage: %prog

Add --help to see a full list of required and optional
arguments to run metaHGT.

Additional information can also be found at:
https://github.com/luo-chengwei/metaHGT/wiki

If you use MeTaxa in your work, please cite it as:
<MeTaxa citation here>

Copyright: Chengwei Luo, Konstantinidis Lab, Georgia Institute of Technology, 2013
"""

import sys
import os
import urllib2


url = "http://enve-omics.ce.gatech.edu/metaxa/db/db.latest.tar.gz"
dir = './'

file_name = dir + url.split('/')[-1]
try:
	u = urllib2.urlopen(url)
except:
	url = url = "http://enve-omics.ce.gatech.edu/mytaxa/db/db.latest.tar.gz"
	u = urllib2.urlopen(url)

f = open(file_name, 'wb')
meta = u.info()
file_size = int(meta.getheaders("Content-Length")[0])
print "Downloading: %s Bytes: %s" % (file_name, file_size)

file_size_dl = 0
block_sz = 8192
while True:
    buffer = u.read(block_sz)
    if not buffer:
        break

    file_size_dl += len(buffer)
    f.write(buffer)
    status = r"%10d  [%3.2f%%]" % (file_size_dl, file_size_dl * 100. / file_size)
    status = status + chr(8)*(len(status)+1)
    print status,
f.close()

os.system('tar xzvf db.latest.tar.gz --directory ../')

