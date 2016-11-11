#!/usr/bin/env python
# --*-- coding:utf-8 --*--
'this is for draw irm graphs together, so that can print directly'

import os
import re
from matplotlib import pyplot as plt
import numpy as np
path = '/backup/jiabo/irmgraph/'

#fig.subplots_adjust(wspace=0.4)
counts = 33

samples = []
index=[]
for line in os.listdir(path):
    if re.search(r'\Spng$', line):
        samples.append(line)
        index.append(line.replace('.irmc.png', '').replace('MSM33-55-1_d', '').replace('_', '.'))

index = [float(i) for i in index]
index_ = np.sort(index)


for i in np.arange(3, 4, step=1):
    fig = plt.figure(figsize=(55, 60), dpi=300, facecolor='white')
    print '--------------------'
    for n in np.arange(6):
        if n+6*i<33:
            print n
            fig.add_subplot(3, 2, n+1)
            fig.subplots_adjust(hspace=0, wspace=0)
            sequen = index_[n+6*i]
            sample = samples[index.index(sequen)]
            im = plt.imread(path+sample)
            plt.imshow(im, interpolation='none')
            plt.axis('off')
            plt.savefig('/home/jiabo/Documents/print/'+'MSM33-55-1_IRMC_'+str(i)+'.png')
