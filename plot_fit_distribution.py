#!/usr/bin/env python3
from numpy.core.fromnumeric import mean
import pandas as pd
import numpy as np
import math
import matplotlib
matplotlib.use("pgf")
import matplotlib.pyplot as plt
plt.rcParams.update({
    "pgf.texsystem": "lualatex",
    "font.family": "serif",  # use serif/main font for text elements
    "text.usetex": True,     # use inline math for ticks
    "pgf.rcfonts": False,    # don't setup fonts from rc parameters
    "axes.labelsize": 20,
    "axes.linewidth": 2.0,
    'axes.unicode_minus': False,
    "font.size": 16,
    "pgf.preamble": '\n'.join([
        "\\usepackage{units}",
        "\\usepackage{metalogo}",
        "\\usepackage{unicode-math}",
        r"\setmathfont{MathJax_Math}",
        r"\setmainfont{FreeSans}",
    ])
})

def kl_divergence(p, q):
    s = 0
    for i in range(len(q)):
        if q[i] == 0.0 or p[i] == 0.0:
            pass
        else:
            s += p[i] * np.log2(p[i] / q[i])
    return s

def PDF_normal(x, mu=0.0, sigma=1.0):
    tmp1 = (x - mu) / sigma
    tmp2 = np.exp(-0.5 * tmp1 * tmp1)
    tmp3 = 1.0 / (sigma * np.sqrt(2.0 * np.pi))
    return tmp3 * tmp2

def PDF_gamma(x, k, theta):
    tmp1 = np.power(x, k - 1.0) * np.exp(-1.0 * x / theta)
    tmp2 = 1.0 / (math.gamma(k) * math.pow(theta, k))
    return tmp1 * tmp2

def PDF_chi(x, k):
    tmp1 = np.power(x, k - 1.0) * np.exp(-1.0 * x * x / 2.0)
    tmp2 = 1.0 / (math.pow(2, 0.5 * k - 1) * math.gamma(0.5 * k))
    return tmp1 * tmp2

def PDF_Maxwell_Boltzmann(x, a):
    tmp1 = x * x * np.exp(-1.0 * x * x / (2.0 * a * a))
    tmp2 = np.sqrt(2.0 / np.pi) / (a * a * a)
    return tmp1 * tmp2

def write_PDF_data(outputfile, x, **kwargs):
    with open(outputfile, 'w') as fOutput:
        fOutput.write('# X')
        for k in kwargs.keys():
            fOutput.write(f' {k}')
        fOutput.write('\n')
        for i in range(len(x)):
            fOutput.write(f'{x[i]:12.7f}')
            for data in kwargs.values():
                fOutput.write(f'{ data[i]:12.7f}')
            fOutput.write('\n')

def recalculate_mean_variance(centers, counts):
    sum_x = 0
    square_sum_x = 0
    count = 0
    for x, y in zip(centers, counts):
        sum_x += x * y
        square_sum_x += x * x * y
        count += y
    mean = sum_x / count
    square_mean = square_sum_x / count
    variance = square_mean - mean * mean
    return mean, variance

data = pd.read_csv('dV_stat.dat', delimiter='\s+', header=None, comment='#')
mean2, var2 = recalculate_mean_variance(data[0], data[1])
print(f'Recalculated mean = {mean2} ; variance = {var2}')
bin_width = data[0][1] - data[0][0]
density = data[1] / sum(data[1]) / bin_width
plt.plot(data[0], density, color='black', label='Data')
# fitting to normal distribution
mu = 4.20394
variance = 20.1803 - mu * mu
print(f'Normal distribution: mu = {mu:.6f} ; sigma = {np.sqrt(variance):.6f}')
normal_curve = PDF_normal(data[0], mu=mu, sigma=np.sqrt(variance))
# fitting to gamma distribution
theta = variance / mu
k = mu / theta
print(f'Gamma distribution: k = {k:.6f} ; theta = {theta:.6f}')
gamma_curve = PDF_gamma(data[0], k=k, theta=theta)
# fitting to chi distribution
k_chi = variance + mu * mu
print(f'Chi distribution: k = {k_chi:.6f}')
chi_curve = PDF_chi(data[0], k_chi)
# fitting to Maxwell-Boltzmann distribution
a = mu / 2.0 / np.sqrt(2.0 / np.pi)
print(f'Maxwell-Boltzmann distribution: a = {a:.6f}')
maxwell_boltzmann_curve = PDF_Maxwell_Boltzmann(data[0], a)
# kl-divergence
kl_pq_normal = kl_divergence(normal_curve, density)
kl_pq_gamma = kl_divergence(gamma_curve, density)
kl_pq_chi = kl_divergence(chi_curve, density)
kl_pq_maxwell_boltzmann = kl_divergence(maxwell_boltzmann_curve, density)
print(f'KL-divergence of fitting by a normal distribution: {kl_pq_normal:.6f}')
print(f'KL-divergence of fitting by a gamma distribution: {kl_pq_gamma:.6f}')
print(f'KL-divergence of fitting by a chi distribution: {kl_pq_chi:.6f}')
print(f'KL-divergence of fitting by a Maxwell-Boltzmann distribution: {kl_pq_maxwell_boltzmann:.6f}')
plt.plot(data[0], normal_curve, label='Normal')
plt.plot(data[0], gamma_curve, label='Gamma')
# plt.plot(data[0], chi_curve, color='green', label='Chi')
plt.plot(data[0], maxwell_boltzmann_curve, label='Maxwell-\nBoltzmann')
plt.legend(prop = {'size': 16}, fancybox = False, frameon = False, handlelength=1.0, handletextpad=0.5)
ax = plt.gca()
ax.xaxis.get_major_formatter()._usetex = False
ax.yaxis.get_major_formatter()._usetex = False
ax.tick_params(direction='in', which='major', length=6.0,
                width=1.0, top=True, right=True, pad=8.0)
ax.tick_params(direction='in', which='minor', length=3.0,
                width=1.0, top=True, right=True, pad=8.0)
ax.set_ylim(0, 0.30)
ax.set_xlim(0, 10)
plt.savefig('dV_stat.png', dpi=300, bbox_inches='tight')
write_PDF_data('density.dat', data[0], density=density,
               p_normal=normal_curve, p_gamma=gamma_curve,
               p_chi=chi_curve, p_maxwell_boltzmann=maxwell_boltzmann_curve)