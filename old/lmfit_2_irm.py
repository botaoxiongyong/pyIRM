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
from scipy.optimize import leastsq
from sklearn.mixture import GMM
from lmfit.models import GaussianModel
from lmfit import minimize, Parameters, Parameter, report_fit
from matplotlib import rcParams

def data():
    """
    将x, gradient(y) 当做是一组密度分布曲线， 根据曲线重新估计随机数列，然后对随机数列进行拟合，并转换坐标
    """

    y_gradient_init = []
    for i in np.gradient(y_measure):
        if i >0:
            y_gradient_init.append(i)
        else:
            y_gradient_init.append(10**-11)

    x_fit = np.linspace(np.log10(x_measure).min(), np.log10(x_measure).max(), 1000)
    y_gradient = interpolate.splev(x_fit, interpolate.splrep(np.log10(x_measure), y_gradient_init))
    data = []
    for i in np.arange(len(x_fit)):
        frequency = y_gradient[i] / sum(y_gradient)
        numbers = np.int(frequency * 10**5)
        if numbers != 0:
            data.extend([x_fit[i]]*numbers)

    data_1D = np.array([[i] for i in np.array(data)])
    return x_fit, y_gradient, data_1D

def loggaussfit():
    x_fit, y_gradient, data_1D = data()
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
        models[i] = GMM(N[i]).fit(data_1D)
    # compute the AIC and the BIC
    AIC = [m.aic(data_1D) for m in models]
    #------------------------------------------------------------
    # Plot the results
    #  We'll use three panels:
    #   1) data + best-fit mixture
    #   2) AIC and BIC vs number of components
    #   3) probability that a point came from each component
    M_best = models[np.argmin(AIC)]
    logprob, responsibilities = M_best.score_samples(x_fit.reshape((-1,1)))#M_best.eval(x)
    pdf = np.exp(logprob)
    pdf_individual = responsibilities * pdf[:, np.newaxis]
    y_fit = (pdf.max())*(y_gradient/(y_gradient.max()))

    return pdf_individual, y_gradient, x_fit, y_fit

def draw_fig():
    params, x_fit, y_gradient, y_fit = lmfit_result()
    fig = plt.figure(figsize=(5, 5), dpi=100, facecolor='white')
    fig.subplots_adjust(left=0.15, right=0.97,
                        bottom=0.15, top=0.9, wspace=0.5, hspace=0.5)
    #fig.suptitle(sample.split('.irmc')[0])
    ax = fig.add_subplot(111, axisbg='white')#'#f7fbff'
    ax.scatter(x_measure, np.gradient(y_measure), facecolors='white', edgecolors='k', s=15, marker='o', alpha=1,
               label='Measured')
    #ax.plot(10**x_fit, y_gradient)
    #ax.plot(10**x_fit, y_fit, color='green')

    y_components = func(x_fit, params)
    fit_curve = sum(y_components)
    p = y_gradient.max()/fit_curve.max()            #转换系数
    x_draw = 10**x_fit
    color = ['#e41a1c','#4daf4a','#377eb8','#984ea3','#ff7f00']
    ax.plot(x_draw, fit_curve*p, '-k', label='Fit')

    label_size = 18
    ax.set_xlabel('Filed (mT)', fontsize=label_size)
    ax.set_ylabel('IRM acquisition gradient', fontsize=label_size)
    ax.set_xscale('log')
    ax.set_xlim(2, 3000)
    ax.set_ylim(0, fit_curve.max()*p*1.1)
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.tick_params(axis='both', which='major', labelsize=label_size)

    y_sum = []
    for i in np.arange(group_number):
        y_sum.append(sum(y_components[i]))

    y_sort =sorted(y_sum, reverse=True)
    for m in np.arange(len(y_sort)):
        for n in np.arange(group_number):
            sequence = 'g'+str(n+1)+'_'
            if y_sort[m] == sum(y_components[n]):
                print m
                x_mean = [x_draw[i] for i in np.arange(len(y_components[n])) if y_components[n][i] == np.max(y_components[n])]
                print x_mean[0], sum(y_components[n])/sum(y_sum)
                content=sum(y_components[n])/sum(y_sum)*100
                label = str(np.int(content))+'%   '+ str(np.int(x_mean[0]))+' mT\n' +\
                        'dp='+str('%.2f'%params[sequence+'sigma'].value)
                ax.plot(x_draw, y_components[n]*p, color=color[m],label=label)
                plt.legend(frameon=False, fontsize = 13)

def lmfit_result():
    pdf_individual, y_gradient, x_fit, y_fit = loggaussfit()
    for i in np.arange(group_number):
        sequence = 'g'+str(i+1)+'_'
        if i ==0:
            gauss = GaussianModel(prefix=sequence)
            pars = gauss.guess(pdf_individual[:,i], x = x_fit)
            center_value = pars[sequence+'center'].value
            pars[sequence+'center'].set(center_value, min=center_value, max=center_value+0.000001)
            sigma_value = pars[sequence+'sigma'].value
            pars[sequence+'sigma'].set(sigma_value, min=sigma_value-0.1, max=sigma_value+0.1)
        else:
            gauss = GaussianModel(prefix=sequence)
            pars.update(gauss.guess(pdf_individual[:,i], x = x_fit))
            center_value = pars[sequence+'center'].value
            pars[sequence+'center'].set(center_value, min=center_value, max=center_value+0.000001)
            sigma_value = pars[sequence+'sigma'].value
            pars[sequence+'sigma'].set(sigma_value, min=sigma_value-0.1, max=sigma_value+0.1)
    params = pars

    result = minimize(residuals, params, args=(x_fit, y_fit), method='cg')
    lmfit_params = result.params
    report_fit(result)
    return lmfit_params, x_fit, y_gradient, y_fit

def lmfit_2():
    lmfit_params, x_fit, y_gradient, y_fit=lmfit_result()
    for i in np.arange(group_number):
        sequence = 'g'+str(i+1)+'_'
        if i ==0:
            gauss = GaussianModel(prefix=sequence)
            pars = gauss.guess(func(x_fit,lmfit_params)[i], x = x_fit)
            center_value = pars[sequence+'center'].value
            pars[sequence+'center'].set(center_value, min=center_value-0.000001, max=center_value+0.000001)
            sigma_value = pars[sequence+'sigma'].value
            pars[sequence+'sigma'].set(sigma_value, min=sigma_value, max=sigma_value+0.00001)
            amplitude_value = lmfit_params[sequence+'amplitude'].value
            pars[sequence+'amplitude'].set(amplitude_value, min=amplitude_value, max=amplitude_value+0.0001)
        else:
            gauss = GaussianModel(prefix=sequence)
            pars.update(gauss.guess(func(x_fit,lmfit_params)[i], x = x_fit))
            center_value = pars[sequence+'center'].value
            pars[sequence+'center'].set(center_value, min=center_value-0.000001, max=center_value+0.000001)
            sigma_value = pars[sequence+'sigma'].value
            pars[sequence+'sigma'].set(sigma_value, min=sigma_value-0.0001, max=sigma_value+0.0001)
            amplitude_value = lmfit_params[sequence+'amplitude'].value
            pars[sequence+'amplitude'].set(amplitude_value, min=amplitude_value, max=amplitude_value+0.0001)
    params = pars

    result = minimize(residuals, params, args=(x_fit, y_fit), method='cg')
    lm2_params = result.params
    report_fit(result)
    return lm2_params, x_fit, y_gradient, y_fit

def leasqr_fit():
    lmfit_params, x_fit, y_gradient, y_fit=lmfit_result()
    params_least=[]
    for i in np.arange(group_number):
        sequence = 'g'+str(i+1)+'_'
        A=lmfit_params[sequence+'amplitude'].value
        s=lmfit_params[sequence+'sigma'].value
        c=lmfit_params[sequence+'center'].value
        params_least.append([A, s, c])
    plsqr = leastsq(residual_least, params_least, args=(x_fit, y_gradient), maxfev=1000)
    plsq = plsqr[0]
    print plsq
    return plsq, x_fit, y_gradient, y_fit

def residual_least(params_least, x, y):
    y_component = func_least(x, params_least)
    return y - sum(y_component)

def func_least(x, params_least):
    y_component=[]
    for i in np.arange(group_number):
        A= params_least[0+3*i]
        s= params_least[1+3*i]
        c= params_least[2+3*i]
        y = A/(s*np.sqrt(2*np.pi))*np.exp(-(x-c)**2/(2*s**2))
        y_component.append(y)
    return y_component

def func(x, params):
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

def residuals(params, x, y):
    y_component = func(x, params)
    return y - sum(y_component)
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
    global y_measure, x_measure, group_number
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
    y_measure = L[1]
    x_measure = L[0]*10**3 #mT

    #lmfit_2()
    draw_fig()
    plt.show()

def save_fig():
    global y_measure, x_measure, group_number
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
            y_measure = L[1]
            x_measure = L[0]*10**3
            draw_fig()
            plt.savefig('/backup/jiabo/irmgraph/'+sample+'.png')

if __name__ == '__main__':
    #main()
    save_fig()
