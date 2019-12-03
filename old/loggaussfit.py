#!/usr/bin/env python
# --*-- coding:UTF-8 --*--
"""
this is a gauss unmixing methos for IRM(isothermal remanent magnetisition) acquisition curves
which is edit from
1D Gaussian Mixture Example
---------------------------
see below and https://github.com/astroML/astroML/blob/master/book_figures/chapter4/fig_GMM_1D.py
"""
import os
import re
from scipy import interpolate
from matplotlib import pyplot as plt
import numpy as np
from sklearn.mixture import GMM


def loggaussfit(x_measure, y):


    """
    将x, gradient(y) 当做是一组密度分布曲线， 根据曲线重新估计随机数列，然后对随机数列进行拟合，并转换坐标
    """

    y_gradient = []
    for i in np.gradient(y):
        if i >0:
            y_gradient.append(i)
        else:
            y_gradient.append(10**-11)

    x_fit = np.log10(x_measure)
    x_interp = np.linspace(x_fit.min(), x_fit.max(), 1000)

    y_interp = interpolate.splev(x_interp, interpolate.splrep(x_fit, y_gradient))

    x_fit = x_interp
    y_gradient = y_interp

    D = []
    for i in np.arange(len(x_fit)):
        xx = x_fit[i]
        yx = y_gradient[i]
        frequency = yx / sum(y_interp)
        numbers = np.int(frequency * 10**5)
        if numbers != 0:
            D.extend([xx]*numbers)

    X = []
    for i in np.array(D):
        X.append([i])
    X = np.array(X)


    #--------------------------------------
    # Learn the best-fit GMM models
    #  Here we'll use GMM in the standard way: the fit() method
    #  uses an Expectation-Maximization approach to find the best
    #  mixture of Gaussians for the data
    # fit models with 1-10 components
    N = np.arange(1, 5)
    models = [None for i in range(len(N))]

    for i in range(len(N)):
        models[i] = GMM(N[i]).fit(X)

    # compute the AIC and the BIC
    AIC = [m.aic(X) for m in models]
    BIC = [m.bic(X) for m in models]

    #------------------------------------------------------------
    # Plot the results
    #  We'll use three panels:
    #   1) data + best-fit mixture
    #   2) AIC and BIC vs number of components
    #   3) probability that a point came from each component

    fig = plt.figure(figsize=(10, 8), dpi=100, facecolor='white')
    fig.subplots_adjust(left=0.12, right=0.97,
                        bottom=0.21, top=0.9, wspace=0.5, hspace=0.5)
    fig.suptitle(sample.split('.irmc')[0])
    color = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00']
    '''
    ----------------------------------------------------------------# plot 1: data + best-fit mixture
    '''
    ax = fig.add_subplot(222)
    M_best = models[np.argmin(AIC)]

    print M_best.params

    x = x_fit#np.linspace(X.min(), X.max(), 1000)#np.array(x_mid)
    logprob, responsibilities = M_best.score_samples(x.reshape((-1,1)))#M_best.eval(x)
    pdf = np.exp(logprob)
    pdf_individual = responsibilities * pdf[:, np.newaxis]

    y_pdf = (np.mean(y_gradient))*(pdf / np.mean(pdf))
    y_in = (np.mean(y_gradient))*(pdf_individual / (np.mean(pdf)))

    #ax.hist(X, 50, normed=True, histtype='stepfilled', alpha=0.4)
    ax.scatter(x_measure, np.gradient(y), facecolors='white', edgecolors='k', s=10, marker='s', alpha=1)
    ax.plot(10**x, y_pdf, '-k')
    for i in np.arange(y_in.shape[1]):
        ax.plot(10**x, y_in[:, i], '--k', color=color[i])
    ax.text(0.04, 0.96, "Best-fit Mixture",
            ha='left', va='top', transform=ax.transAxes)
    ax.set_xlabel('Filed (mT)')
    ax.set_ylabel('IRM acquisition gradient')
    ax.set_ylim(0, 10**-8)
    ax.set_xscale('log')
    ax.set_xlim(2, 3000)
    ax.set_ylim(-y_pdf.max()*0.1, y_pdf.max()*1.3)
    '''
    ---------------------------------------------------------------##plot add: irm aqcuisition curve:
    '''
    ax = fig.add_subplot(221)
    ax.scatter(x_measure, y/y.max(), facecolors='white', edgecolors='k', s=10, marker='s', alpha=0.5)
    y_irm = y_pdf / sum(y_pdf)
    y_acumulation = []
    m = 0
    for i in y_irm:
        m = m+i
        y_acumulation.append(m)
    ax.plot(10**x, y_acumulation, '-k')

    y_irm = y_in / sum(y_pdf)
    for i in np.arange(y_irm.shape[1]):
        y_acumulation = []
        m = 0
        for n in y_irm[:, i]:
            m = m + n
            y_acumulation.append(m)
        ax.plot(10**x, y_acumulation, '--k', color=color[i])
    ax.set_xlabel('Filed (mT)')
    ax.set_ylabel('normalized IRM')
    ax.set_xscale('log')
    ax.set_xlim(2, 3000)
    ax.set_ylim(0, 1.1)
    '''
    --------------------------------------------------------------# plot 2: AIC and BIC
    '''
    ax = fig.add_subplot(223)
    ax.plot(N, AIC, '-r', label='AIC', alpha=0.5)
    ax.plot(N, BIC, '--k', label='BIC')
    ax.set_xlabel('n. components')
    ax.set_ylabel('information criterion')
    ax.legend(loc=1, frameon=False)
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    '''
    --------------------------------------------------------------# plot 3: posterior probabilities for each component
    '''
    ax = fig.add_subplot(224)
    p = M_best.predict_proba(x.reshape((-1,1)))#M_best.predict_proba(x)
    #p = p[:, (1, 0, 2)]  # rearrange order so the plot looks better
    p = p.cumsum(1).T
    ax.fill_between(10**x, 0, p[0], color=color[0], alpha=0.3)
    for i in np.arange(1, p.shape[0]):
        if i < p.shape[0] - 1:
            ax.fill_between(10**x, p[i-1], p[i], color=color[i], alpha=0.3)
        if i == p.shape[0] - 1:
            ax.fill_between(10**x, p[i-1], 1, color=color[i], alpha=0.3)
    #ax.set_xlim(0, 4)
    ax.set_ylim(0, 1)
    ax.set_xlabel('Filed (mT)')
    ax.set_ylabel(r'posterior probabilities')
    ax.set_xscale('log')
    ax.set_xlim(3, 2000)


    #ax.text(0, 3, 'class 1', rotation='vertical')
    #ax.text(0.3, 3, 'class 2', rotation='vertical')
    #ax.text(3, 5, 'class 3', rotation='vertical')


path = '/backup/jiabo/rockdata/MSM33-55-1/'
sample = 'MSM33-55-1_d99.irmc'
l = []
with open(path + sample) as dat:
    data = dat.readlines()
    data = data[87:-2]
    for lines in data:
        var = lines.split(',')
        l.append(var)
L = np.array(l, dtype=np.float64).T
y = L[1]
x_measure = L[0]*10**3

loggaussfit(x_measure, y)

plt.show()

'''
for line in os.listdir(path):
    if re.search(r'\Sirmc$', line):
        sample = line.replace('.irmc', '')
        loggaussfit(path, sample+'.irmc')
        plt.savefig('/backup/jiabo/irmgraph/'+sample)
'''

'''
"""
1D Gaussian Mixture Example
---------------------------
Figure 4.2.
Example of a one-dimensional Gaussian mixture model with three components.
The left panel shows a histogram of the data, along with the best-fit model
for a mixture with three components. The center panel shows the model selection
criteria AIC (see Section 4.3) and BIC (see Section 5.4) as a function of the
number of components. Both are minimized for a three-component model. The
right panel shows the probability that a given point is drawn from each class
as a function of its position. For a given x value, the vertical extent of
each region is proportional to that probability. Note that extreme values
are most likely to belong to class 1.
"""
# Author: Jake VanderPlas
# License: BSD
#   The figure produced by this code is published in the textbook
#   "Statistics, Data Mining, and Machine Learning in Astronomy" (2013)
#   For more information, see http://astroML.github.com
#   To report a bug or issue, use the following forum:
#    https://groups.google.com/forum/#!forum/astroml-general
from matplotlib import pyplot as plt
import numpy as np
from sklearn.mixture import GMM

#----------------------------------------------------------------------
# This function adjusts matplotlib settings for a uniform feel in the textbook.
# Note that with usetex=True, fonts are rendered with LaTeX.  This may
# result in an error if LaTeX is not installed on your system.  In that case,
# you can set usetex to False.
from astroML.plotting import setup_text_plots
setup_text_plots(fontsize=8, usetex=True)

#------------------------------------------------------------
# Set up the dataset.
#  We'll use scikit-learn's Gaussian Mixture Model to sample
#  data from a mixture of Gaussians.  The usual way of using
#  this involves fitting the mixture to data: we'll see that
#  below.  Here we'll set the internal means, covariances,
#  and weights by-hand.
np.random.seed(1)

gmm = GMM(3, n_iter=1)
gmm.means_ = np.array([[-1], [0], [3]])
gmm.covars_ = np.array([[1.5], [1], [0.5]]) ** 2
gmm.weights_ = np.array([0.3, 0.5, 0.2])

X = gmm.sample(1000)

#------------------------------------------------------------
# Learn the best-fit GMM models
#  Here we'll use GMM in the standard way: the fit() method
#  uses an Expectation-Maximization approach to find the best
#  mixture of Gaussians for the data

# fit models with 1-10 components
N = np.arange(1, 11)
models = [None for i in range(len(N))]

for i in range(len(N)):
    models[i] = GMM(N[i]).fit(X)

# compute the AIC and the BIC
AIC = [m.aic(X) for m in models]
BIC = [m.bic(X) for m in models]

#------------------------------------------------------------
# Plot the results
#  We'll use three panels:
#   1) data + best-fit mixture
#   2) AIC and BIC vs number of components
#   3) probability that a point came from each component

fig = plt.figure(figsize=(5, 1.7))
fig.subplots_adjust(left=0.12, right=0.97,
                    bottom=0.21, top=0.9, wspace=0.5)


# plot 1: data + best-fit mixture
ax = fig.add_subplot(131)
M_best = models[np.argmin(AIC)]

x = np.linspace(-6, 6, 1000)
logprob, responsibilities = M_best..score_samples(x.reshape((-1,1)))#eval(x)
pdf = np.exp(logprob)
pdf_individual = responsibilities * pdf[:, np.newaxis]

ax.hist(X, 30, normed=True, histtype='stepfilled', alpha=0.4)
ax.plot(x, pdf, '-k')
ax.plot(x, pdf_individual, '--k')
ax.text(0.04, 0.96, "Best-fit Mixture",
        ha='left', va='top', transform=ax.transAxes)
ax.set_xlabel('$x$')
ax.set_ylabel('$p(x)$')


# plot 2: AIC and BIC
ax = fig.add_subplot(132)
ax.plot(N, AIC, '-k', label='AIC')
ax.plot(N, BIC, '--k', label='BIC')
ax.set_xlabel('n. components')
ax.set_ylabel('information criterion')
ax.legend(loc=2)


# plot 3: posterior probabilities for each component
ax = fig.add_subplot(133)

p = M_best.predict_proba(x.reshape((-1,1)))#predict_proba(x)
p = p[:, (1, 0, 2)]  # rearrange order so the plot looks better
p = p.cumsum(1).T

ax.fill_between(x, 0, p[0], color='gray', alpha=0.3)
ax.fill_between(x, p[0], p[1], color='gray', alpha=0.5)
ax.fill_between(x, p[1], 1, color='gray', alpha=0.7)
ax.set_xlim(-6, 6)
ax.set_ylim(0, 1)
ax.set_xlabel('$x$')
ax.set_ylabel(r'$p({\rm class}|x)$')

ax.text(-5, 0.3, 'class 1', rotation='vertical')
ax.text(0, 0.5, 'class 2', rotation='vertical')
ax.text(3, 0.3, 'class 3', rotation='vertical')

plt.show()
'''
