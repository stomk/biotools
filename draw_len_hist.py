#!/usr/bin/env python

# input: output file of 'fatt len'

import sys
import os.path
import matplotlib as mpl
import matplotlib.pyplot as plt
from pylab import axes

# input
fin_path = os.path.abspath(sys.argv[1])
fin = open(fin_path, 'r')

vals = [float(line.strip('\n')) for line in fin]
num_reads = len(vals)

fin.close()

# output
name = os.path.splitext(fin_path)[0]
fout_path = name + '.len-hist.png'

# figure
fig = plt.figure()
ax = fig.add_subplot(111)
axes(axisbg='#eeeeee')

#plt.locator_params(axis='x', nbins=10)
plt.hist(vals, range=(0,50000), bins=25, rwidth=0.9, facecolor='g')

plt.xlabel('Length')
plt.ylabel('# reads')
plt.title(name)

num_str = 'total reads: ' + ('%d' % num_reads)
plt.text(0.76, 0.95, num_str, fontsize=15, transform=ax.transAxes)

plt.savefig(fout_path)
