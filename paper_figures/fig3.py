import copy
import glob
import itertools
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerTuple
import numpy as np
import os
import scipy.constants as c
from scipy.interpolate import interp1d
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
from scipy.special import voigt_profile
import sys
# from libs import LindbladMasterEq as ats
import LindbladMasterEq as ats
from elecsus.libs import BasisChanger as bc

os.chdir('paper_figures')
E_LCP = bc.lrz_to_xyz(np.array([1,0,0]))

def fit(x, x0, a, b):
    return a * voigt_profile(x+x0, b, 6.0) + 1


if __name__ == '__main__':
    propagation = False
    # transit_type='integral'
    transit_type='single'
    doppler = True
    x = np.linspace(5500, 18000, 2000)

    # diameter, temperature, buffer gas, intensity correction, unused variable
    pairs=[
        (0.5e-3, 62.0, 0, 0.25, 1),
        (1.0e-3, 66.5, 0, 0.83, 1),
        (2.0e-3, 67.0, 0, 0.90, 1),
        ]

    factors = [
        [0.274, 8.57, 152.5], # 0.5 mm
        [1.270, 12.7, 228.0], # 1.0 mm
        [0.343, 11.0, 343.4], # 2.0 mm
    ]

    # Uncommenting line 131 allows to correct for small temperature fluctuations
    # 1.1 corresponds to 10% which are around 1 degree
    scale_factors = [
        [1.00, 1.01, 1.75 ],
        [1.05, 1.00, 1.00 ],
        [1.09, 1.00, 1.10 ],
    ]


    folders = ['data/0.5mm_0.6T_15V',
               'data/1.0mm_0.6T_15V',
               'data/2.0mm_0.6T_15V',
    ]
    folder_prefix = 'poly7'
    files = [
        ['1.8uW.npz', '56uW.npz', '1mW.npz'], # 0.5 mm
        ['10uW.npz', '100uW.npz', '1.8mW.npz'], # 1.0 mm
        ['10uW.npz', '320uW.npz', '10mW.npz'], # 2.0 mm
    ]

    cmap = plt.get_cmap('viridis', 5)
    def cmap(i):
        return f'C{i}'
    plt.rcParams['axes.prop_cycle'] = plt.cycler(color=plt.cm.tab20.colors)
    lot1 = mpl.lines.Line2D([], [], ls='-', c='C0', mfc='none', mec='k', label=r'$\lessapprox$ 1$\cdot$I$_{sat}$')
    lot2 = mpl.lines.Line2D([], [], ls='-', c='C1', mfc='none', mec='k', label=r'$\approx$10$\cdot$I$_{sat}$')
    lot3 = mpl.lines.Line2D([], [], ls='-', c='C2', mfc='none', mec='k', label=r'$\gtrapprox$100$\cdot$I$_{sat}$')
    lct = mpl.lines.Line2D([], [], ls='-', c='k', label='simulation')

    fig0, ax = plt.subplots(2, 3, sharex=True, sharey='row', gridspec_kw={'height_ratios': [3, 1]}, figsize=(8,4))

    for i, pair in enumerate(pairs):
        D = pair[0]
        T = pair[1]
        Isat = 16.69 # W/m2
        Psat = np.pi * Isat * D**2 / 8

        fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1]})
        for j, factor in enumerate(factors[i]):
            print(f'Pair: {pair}; corresponds to data at {Psat * factor / pair[3] * 1e6} uW')
            with np.load(f'{folders[i]}/{folder_prefix}/{files[i][j]}') as f:
                x_exp = f['x']
                y_exp = f['y']
                # Find peak of left and right most peaks to scale to simulation
                ytmp = (y_exp - y_exp.min()) / y_exp.ptp()
                idx_peaks, _ = find_peaks(1-ytmp, 1-0.92, distance=2000 * np.diff(x_exp)[0], prominence=0.05)
                xC = np.array([x_exp[idx_peaks[0]] - 200, x_exp[idx_peaks[0]] + 200])
                xO = np.array([x_exp[idx_peaks[-1]] - 200, x_exp[idx_peaks[-1]] + 200])
                pc, _ = curve_fit(fit, x_exp[(x_exp>xC[0]) & (x_exp<xC[1])], y_exp[(x_exp>xC[0]) & (x_exp<xC[1])], p0=(-x_exp[idx_peaks[0]], -200, 280))
                po, _ = curve_fit(fit, x_exp[(x_exp>xO[0]) & (x_exp<xO[1])], y_exp[(x_exp>xO[0]) & (x_exp<xO[1])], p0=(-x_exp[idx_peaks[-1]], -100, 280))
                xminC = x_exp[np.argmin(fit(x_exp, *pc))]
                xminO = x_exp[np.argmin(fit(x_exp, *po))]
                x_distance_exp = xminO - xminC
                x_distance_theory = 17112 - 6087.79
                x_scale = x_distance_exp / x_distance_theory
                x0 = xminC / x_scale - 6087.79
                x_exp = x_exp / x_scale
                x_exp = x_exp - x0
                exp_interp = interp1d(x_exp, y_exp, kind='cubic', fill_value=1, )
                y_exp = exp_interp(x)

            lcell = 2e-3
            file = f'tmp_fig3/{pair}_{factor}.npz'
            if not os.path.exists(file):
                p_dict = {'Elem':'Rb','Dline':'D2', 'lcell': lcell, 'T': pair[1],
                        'Constrain': True, 'Bfield': 6000, 'rb85frac': 0,
                        'laserPower': 1e-15, 'laserWaist': pair[0],
                        'GammaBuf': pair[2], 'collisions': 'dephase'}
                atom = ats.atomicSystem('Rb87', E_in=E_LCP, p_dict=p_dict)
                T = atom.transmission([ats.beam(w=x, P=factor * Psat, D=D, profile='flat')],
                                    z=lcell, doppler=doppler, transit_type=transit_type)
                np.savez(file, x=x, y=T)
            else:
                with np.load(file) as f:
                    T = f['y']

            # Uncomment to correct for small temperature fluctuations
            # T = np.exp(scale_factors[i][j] * np.log(T))

            # Position of all 8 HPB resonance peaks
            # Used for better linearisation of the frequency axis
            n = [ 6087.79389695,  7238.36918459,  8582.7913957,  10252.37618809, 11496.74837419, 13841.67083542, 15630.06503252, 17112.05602801]
            x_exp_peaks = np.zeros(8)
            x_sim_peaks = np.zeros(8)
            for ni, nval in enumerate(n):
                xfit = x[(x > nval-300) & (x < nval+300)]
                yfitE = y_exp[(x > nval-300) & (x < nval+300)]
                yfitS = T[(x > nval-300) & (x < nval+300)]
                poptE, _ = curve_fit(fit, xfit, yfitE, p0=(-nval, -200, 280))
                poptS, _ = curve_fit(fit, xfit, yfitS, p0=(-nval, -200, 280))
                x_exp_peaks[ni] = x[np.argmin(fit(x, *poptE))]
                x_sim_peaks[ni] = x[np.argmin(fit(x, *poptS))]
            linearized = np.poly1d(np.polyfit(x_exp_peaks, x_sim_peaks, deg=6))
            xnew = np.linspace(linearized(x).min(), linearized(x).max(), 2000)
            y_exp = interp1d(linearized(x), y_exp, kind='cubic', fill_value='extrapolate')(xnew)
            T = interp1d(x, T, kind='cubic', fill_value='extrapolate')(xnew)

            xnew /= 1e3 # make it GHz
            ax1.plot(xnew, y_exp, c=cmap(j), label=factor)
            ax1.plot(xnew, T, '--', c=cmap(j))
            ax2.plot(xnew, (T-y_exp)*100, c=cmap(j))
            ax[0, i].plot(xnew, T, c='k')
            ax[0, i].plot(xnew, y_exp, '-', c=cmap(j), label=factor)
            ax[1, i].plot(xnew, (T-y_exp)*100, c=cmap(j))
        ax1.set_ylabel('Transmission')
        ax2.set_ylabel('Residual [%]')
        ax2.set_xlabel('Detuning [GHz]')
        ax[1, i].set_xlabel('Detuning [GHz]')
        ax[1, i].set_xlim(xnew.min(), xnew.max())
        h, _ = ax1.get_legend_handles_labels()
        fig.tight_layout()
        fig.subplots_adjust(hspace=0, wspace=0)
        plt.savefig(f'tmp_fig3/{pair}_paper.png', dpi=200)
        plt.close()
    ax[0, 0].set_ylabel('Transmission')
    ax[1, 0].set_ylabel('Residual [%]')
    ax[0, 0].text(12000, 0.55, '$⌀=0.5$ mm')
    ax[0, 1].text(12000, 0.55, '$⌀=1.0$ mm')
    ax[0, 2].text(12000, 0.55, '$⌀=2.0$ mm')
    fig0.legend(ncol=4, handles=[lot1, lot2, lot3, lct], frameon=False, loc='upper center', bbox_to_anchor=(0.5, 1.02))
    fig0.subplots_adjust(hspace=0, wspace=0, top=0.93)
    fig0.tight_layout(pad=0)
    fig0.savefig(f'fig3.png', dpi=200)
    plt.show()
