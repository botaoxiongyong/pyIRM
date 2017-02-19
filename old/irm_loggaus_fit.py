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
from lmfit.models import GaussianModel
from lmfit import minimize, Parameters, Parameter, report_fit

def data(x_measure, y):
    """
    将x, gradient(y) 当做是一组密度分布曲线， 根据曲线重新估计随机数列，然后对随机数列进行拟合，并转换坐标
    """

    y_gradient_init = []
    for i in np.gradient(y):
        if i >0:
            y_gradient_init.append(i)
        else:
            y_gradient_init.append(10**-11)

    x_ = np.log10(x_measure)
    x_fit = np.linspace(x_.min(), x_.max(), 1000)
    y_gradient = interpolate.splev(x_fit, interpolate.splrep(x_, y_gradient_init))
    D = []
    for i in np.arange(len(x_fit)):
        xx = x_fit[i]
        yx = y_gradient[i]
        frequency = yx / sum(y_gradient)
        numbers = np.int(frequency * 10**5)
        if numbers != 0:
            D.extend([xx]*numbers)
    X = []
    for i in np.array(D):
        X.append([i])
    X = np.array(X)
    return x_fit, y_gradient, X


def loggaussfit(sample, x_measure, y, group_number):
    x_fit, y_gradient, X = data(x_measure, y)

    #--------------------------------------
    # Learn the best-fit GMM models
    #  Here we'll use GMM in the standard way: the fit() method
    #  uses an Expectation-Maximization approach to find the best
    #  mixture of Gaussians for the data
    # fit models with 1-10 components

    #N = np.arange(group_number, group_number+1)
    N = np.arange(group_number, group_number+1)
    models = [None for i in range(len(N))]

    for i in range(len(N)):
        models[i] = GMM(N[i]).fit(X)

    # compute the AIC and the BIC
    AIC = [m.aic(X) for m in models]
    #------------------------------------------------------------
    # Plot the results
    #  We'll use three panels:
    #   1) data + best-fit mixture
    #   2) AIC and BIC vs number of components
    #   3) probability that a point came from each component



    M_best = models[np.argmin(AIC)]
    x = x_fit#np.linspace(X.min(), X.max(), 1000)#np.array(x_mid)
    logprob, responsibilities = M_best.score_samples(x.reshape((-1,1)))#M_best.eval(x)
    pdf = np.exp(logprob)
    pdf_individual = responsibilities * pdf[:, np.newaxis]
    y_fit = (pdf.max())*(y_gradient/(y_gradient.max()))

    fig = plt.figure(figsize=(5, 5), dpi=100, facecolor='white')
    fig.subplots_adjust(left=0.12, right=0.97,
                        bottom=0.21, top=0.9, wspace=0.5, hspace=0.5)
    fig.suptitle(sample.split('.irmc')[0])
    ax = fig.add_subplot(111, axisbg='white')#'#f7fbff'
    ax.scatter(x_measure, np.gradient(y), facecolors='white', edgecolors='k', s=8, marker='o', alpha=1, label='Measured')

    lmfit_result(pdf_individual, y_fit, ax, x, group_number, y)

def lmfit_result(pdf_individual, y_fit, ax, x, group_number, y):
    for i in np.arange(group_number):
        sequence = 'g'+str(i+1)+'_'
        if i ==0:
            gauss = GaussianModel(prefix=sequence)
            pars = gauss.guess(pdf_individual[:,i], x = x)
            center_value = pars[sequence+'center'].value
            pars[sequence+'center'].set(center_value, min=center_value, max=center_value+0.00001)
            sigma_value = pars[sequence+'sigma'].value
            pars[sequence+'sigma'].set(sigma_value, min=sigma_value-0.1, max=sigma_value+0.1)
        else:
            gauss = GaussianModel(prefix=sequence)
            pars.update(gauss.guess(pdf_individual[:,i], x = x))
            center_value = pars[sequence+'center'].value
            pars[sequence+'center'].set(center_value, min=center_value, max=center_value+0.00001)
            sigma_value = pars[sequence+'sigma'].value
            pars[sequence+'sigma'].set(sigma_value, min=sigma_value-0.1, max=sigma_value+0.1)
    params = pars

    result = minimize(residuals, params, args=(x, y_fit, group_number), method='cg')
    report_fit(result.params)

    y_components = func(x, result.params, group_number)
    fit_curve = sum(y_components)
    p = np.gradient(y).max()/fit_curve.max()            #转换系数

    x_draw = 10**x
    color = ['#e41a1c','#4daf4a','#377eb8','#984ea3','#ff7f00']
    ax.plot(x_draw, fit_curve*p, '-k', label='Fit')
    ax.set_xlabel('Filed (mT)')
    ax.set_ylabel('IRM acquisition gradient')
    ax.set_xscale('log')
    ax.set_xlim(2, 3000)
    ax.set_ylim(-0.1*p, fit_curve.max()*p*1.1)
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

    y_sum = []
    for i in np.arange(group_number):
        y_sum.append(sum(y_components[i]))

    y_sort =sorted(y_sum, reverse=True)
    for m in np.arange(len(y_sort)):
        for n in np.arange(group_number):
            if y_sort[m] == sum(y_components[n]):
                print m
                x_mean = [x_draw[i] for i in np.arange(len(y_components[n])) if y_components[n][i] == np.max(y_components[n])]
                print x_mean[0], sum(y_components[n])/sum(y_sum)
                content=sum(y_components[n])/sum(y_sum)*100
                label = str(np.int(content))+'%   '+ str(np.int(x_mean[0]))+' mT'
                ax.plot(x_draw, y_components[n]*p, color=color[m],label=label)
                plt.legend(frameon=False, fontsize = 11)




def func(x, params, group_number):
    #p, p1, p2, a, a1, a2, b, b1, b2 = para\
    y_component=[]
    for i in np.arange(group_number):
        A = 'g'+str(i+1)+'_amplitude'
        s = 'g'+str(i+1)+'_sigma'
        c = 'g'+str(i+1)+'_center'
        A = params[A].value
        s= params[s].value
        c= params[c].value
        y = A/(s*np.sqrt(2*np.pi))*np.exp(-(x-c)**2/(2*s**2))
        y_component.append(y)
    return y_component

def residuals(params, x, y_fit, group_number):
    y_component = func(x, params, group_number)
    return y_fit - sum(y_component)

#
#single is for single endmember result
def single(x, params, i):
    A = 'g'+str(i+1)+'_amplitude'
    s = 'g'+str(i+1)+'_sigma'
    c = 'g'+str(i+1)+'_center'
    A = params[A].value
    s= params[s].value
    c= params[c].value
    y = A/(s*np.sqrt(2*np.pi))*np.exp(-(x-c)**2/(2*s**2))
    return y



def main():
    path = '/backup/jiabo/rockdata/MSM33-55-1/'
    sample = 'MSM33-55-1_d80_2.irmc'
    group_number = 3
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
    loggaussfit(sample, x_measure, y, group_number)
    plt.show()

def save_fig():
    path = '/backup/jiabo/rockdata/MSM33-55-1/'
    group_number = 3
    for line in os.listdir(path):
        l = []
        if re.search(r'\Sirmc$', line):
            sample = line
            print sample
            with open(path + sample) as dat:
                data = dat.readlines()
                data = data[87:-2]
                for lines in data:
                    var = lines.split(',')
                    l.append(var)
            L = np.array(l, dtype=np.float64).T
            y = L[1]
            x_measure = L[0]*10**3
            loggaussfit(sample, x_measure, y, group_number)

            plt.savefig('/backup/jiabo/irmgraph/'+sample+'.png')

if __name__ == '__main__':
    main()
    #save_fig()
