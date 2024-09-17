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
E_RCP = bc.lrz_to_xyz(np.array([0,1,0]))

groundState = ats.state(5, 0, 1/2)     # 5S1/2
excitedState_D2 = ats.state(5, 1, 3/2)  # 5P1/2
basefolder = 'data'
diameters = np.array([0.5, 1, 2])
# cmap = plt.get_cmap('flare', diameters.size+1)
mss, msl, sym = 8, 8, '.'
offcolor = 'darkgrey'


def make_inv(fig):
    for ax in fig.get_axes():
        ax.tick_params(axis='both', colors='w', which='both')
        ax.title.set_color('w')
        legend = ax.get_legend()
        if legend:
            plt.setp(legend.get_texts(), color='w')
        ax.xaxis.label.set_color('w')
        ax.yaxis.label.set_color('w')
        for spine in ax.spines.values():
            spine.set_edgecolor('w')

# https://stackoverflow.com/a/71887460/7836757
def add_headers(
    fig,
    *,
    row_headers=None,
    col_headers=None,
    row_pad=1,
    col_pad=5,
    rotate_row_headers=True,
    **text_kwargs
):
    # Based on https://stackoverflow.com/a/25814386

    axes = fig.get_axes()

    for ax in axes:
        sbs = ax.get_subplotspec()

        # Putting headers on cols
        if (col_headers is not None) and sbs.is_first_row():
            ax.annotate(
                col_headers[sbs.colspan.start],
                xy=(0.5, 1),
                xytext=(0, col_pad),
                xycoords="axes fraction",
                textcoords="offset points",
                ha="center",
                va="baseline",
                **text_kwargs,
            )

        # Putting headers on rows
        if (row_headers is not None) and sbs.is_first_col():
            ax.annotate(
                row_headers[sbs.rowspan.start],
                xy=(0, 0.5),
                xytext=(-ax.yaxis.labelpad - row_pad, 0),
                xycoords=ax.yaxis.label,
                textcoords="offset points",
                ha="right",
                va="center",
                rotation=rotate_row_headers * 90,
                **text_kwargs,
            )


def fit(x, x0, a, b):
    return a * voigt_profile(x+x0, b, 6.0) + 1


def increase_ijk():
    # increases first k, then j, then i
    global i, j, k
    k += 1
    if k >= k_dim:
        k = 0
        j += 1
        if j >= j_dim:
            j = 0
            i += 1

def get_exp_peaks(set, bfield=False):
    files = glob.glob(f'{set}/*.npz')
    powers = np.empty((len(files)))
    for i, f in enumerate(files):
        powers[i] = float(f.split(f'/')[-1].split('W')[0][:-1])
        if f.split(f'/')[-1].split('W')[0][-1] == 'u':
            powers[i] *= 1e-6
        else:
            powers[i] *= 1e-3
    powers *= 1e6 # uW
    D = float(files[0].split('/')[1].split('mm')[0]) / 10 # cm

    idx = np.argsort(powers)
    powers = powers[idx]
    Isat = 1.669e3 # uW / cm2
    P2I = lambda P: 2 * P / np.pi / (D/2)**2
    I = P2I(powers) / Isat
    lcell = 2e-3

    if bfield == False:
        TF1 = np.zeros_like(powers)
        TF2 = np.zeros_like(powers)
        for i, file in enumerate(np.asarray(files)[idx]):
            with np.load(file) as f:
                x = f['x']
                y = f['y']
            x = x - x[np.argmin(y)] - 2500
            poptF1, _, = curve_fit(fit, x[x>0], y[x>0], p0=(-4000, -350, 278))
            poptF2, _, = curve_fit(fit, x[x<0], y[x<0], p0=(+2500, -350, 278))
            TF1[i] = np.min(fit(x, *poptF1))
            TF2[i] = np.min(fit(x, *poptF2))
        alphaF1 = -np.log(TF1) / lcell
        alphaF2 = -np.log(TF2) / lcell
        alphaF1 /= alphaF1.max()
        alphaF2 /= alphaF2.max()
        return I, alphaF1, alphaF2
    else:
        TC = np.empty(powers.size)
        TO = np.empty(powers.size)
        bc = np.array([(0, 2250), (2250, 3500), (3500, 5100), (5100, 6600)])
        bo = np.array([(6600, 8200), (8200, 10300), (10300, 12000), (12000, 14000)])
        for i, file in enumerate(np.asarray(files)[idx]):
            with np.load(file) as f:
                x = f['x']
                y = f['y']
            # We shift the spectrum to find the left most peak.
            # However, it is not always the strongest, therefore we search at higher transmissions
            ptp = 1 - y.min()
            arg = 1 - 0.9 * ptp
            x = x - x[np.where(y<arg)[0][0]] + 1500
            try:
                argminC = np.argmin(y[(x>bc[0][0]) & (x<bc[0][1])]) + np.argwhere(x>bc[0][0]).min()
                argminO = np.argmin(y[(x>bo[3][0]) & (x<bo[3][1])]) + np.argwhere(x>bo[3][0]).min()
                xC = np.array([x[argminC] - 100, x[argminC] + 100])
                xO = np.array([x[argminO] - 100, x[argminO] + 100])
                pc, _ = curve_fit(fit, x[(x>xC[0]) & (x<xC[1])], y[(x>xC[0]) & (x<xC[1])], p0=(-x[argminC], -700*ptp, 280))
                po, _ = curve_fit(fit, x[(x>xO[0]) & (x<xO[1])], y[(x>xO[0]) & (x<xO[1])], p0=(-x[argminO], -275*ptp, 280))
                TC[i] = np.min(fit(x, *pc))
                TO[i] = np.min(fit(x, *po))
            except:
                print(f'Failed for {file} in set {set}')
                TC[i] = np.nan
                TO[i] = np.nan
        alphaC = -np.log(TC) / lcell
        alphaO = -np.log(TO) / lcell
        if normalized:
            alphaC /= alphaC.max()
            alphaO /= alphaO.max()
        return I, alphaC, alphaO

def get_sim_peaks(folder, pair, suffix):
    filename = f'{folder}{pair[0]*1000:.1f}_{pair[1]:.1f}_{pair[2]:.1f}{suffix}.npz'
    if os.path.exists(filename):
        with np.load(filename) as f:
            Isat_factor = f['Isat_factor']
            alphaO = f['alphaO']
            alphaC = f['alphaC']
    else:
        print(f'Warning, could not find file for pair {folder}{pair[0]*1000:.1f}_{pair[1]:.1f}_{pair[2]:.1f}{suffix}.npz')
        Isat_factor = np.geomspace(1e-3, 1e3, 15)
        alphaO = np.zeros_like(Isat_factor)
        alphaC = np.zeros_like(Isat_factor)
    if normalized:
        alphaC /= alphaC.max()
        alphaO /= alphaO.max()
    return Isat_factor, alphaC, alphaO


if __name__ == '__main__':
    # #########################################################################
    # # 0.6 T
    # #########################################################################
    # Construct 3d array of folder names
    # with dimensions: [temperature, diameter, bufgas]
    # we however also have the open/closed transition differenciation
    # the 3/5/7 denominator indicates the highest polynomial order used for frequency correction and transmission normalisation
    filenames3 = np.empty((3, 2, 2), dtype="object")
    filenames5 = np.empty((3, 2, 2), dtype="object")
    filenames7 = np.empty((3, 2, 2), dtype="object")
    for di, d in enumerate(['0.5mm', '1.0mm', '2.0mm']):
        for ti, t in enumerate(['15V', '17.5V']):
            for gi, g in enumerate(['', '_bufgas']):
                filenames3[di,ti,gi] = f'{basefolder}/{d}_0.6T_{t}{g}/poly3/'
                filenames5[di,ti,gi] = f'{basefolder}/{d}_0.6T_{t}{g}/poly5/'
                filenames7[di,ti,gi] = f'{basefolder}/{d}_0.6T_{t}{g}/poly7/'
    # For 6000G, the peak absorption happens at 6087,79 and 17112.05 Mhz. It does not significantly change,
    # therefore we can track absorption at those single values during propagation
    N = 15
    i_dim = 3
    j_dim = 2
    k_dim = 1
    pairs = np.array([
        (0.5e-3, 62.0, 0, 0.25, 1),
        (0.5e-3, 81.0, 0, 0.28, 1),
        (1.0e-3, 66.5, 0, 0.83, 1),
        (1.0e-3, 81.5, 0, 0.88, 1.01),
        (2.0e-3, 67.0, 0, 0.90, 1),
        (2.0e-3, 81.5, 0, 0.85, 1.02),
    ]).reshape((i_dim, j_dim, k_dim, -1))
    normalized = False

    # Some shared variables
    labelT = ['low T', 'high T']
    labelG = ['w/o buffergas', 'w/ buffergas']
    labelD = ['diameter = 0.5 mm', 'diameter = 1.0 mm ', 'diameter = 2.0 mm']
    labelD_alt = [r'$⌀=0.5$ mm', r'$⌀=1.0$ mm', r'$⌀=2.0$ mm']
    lot = mpl.lines.Line2D([], [], marker='o', ls='--', c='k', mfc='none', mec='k', label='open transition')
    lct = mpl.lines.Line2D([], [], marker='o', ls='-', c='k', label='closed transition')
    font_kwargs = dict(fontfamily="serif", fontsize="large")

    cmap = plt.get_cmap('coolwarm', 2)
    fig, ax = plt.subplots(1, 3, sharex=True, sharey=True, figsize=(10, 3.5))
    fig.subplots_adjust(hspace=0, wspace=0)
    i = 0  # beam diameter
    j = 0  # temperature
    k = 0  # buffer gas
    for p, pair in enumerate((pairs.reshape(-1, 5))):
        data = get_exp_peaks(filenames7[i, j, k], bfield=True)
        ax[i].scatter(data[0]*pair[3], data[1]*pair[4], c=cmap(j), label=labelT[j])
        ax[i].scatter(data[0]*pair[3], data[2]*pair[4], facecolors='none', edgecolors=cmap(j))
        sim = get_sim_peaks('simulation/', pair, '_transit_single_factor_propagated')
        ax[i].plot(sim[0], sim[1], c=cmap(j), ls='-')
        ax[i].plot(sim[0], sim[2], c=cmap(j), ls='--')
        increase_ijk()

    [ax.set_xscale('log') for ax in fig.get_axes()]
    [ax[i].text(20, 900, labelD_alt[i], fontsize=11) for i in range(3)]
    ax[0].set_xlabel(r'$I/I_{sat}$')
    ax[1].set_xlabel(r'$I/I_{sat}$')
    ax[2].set_xlabel(r'$I/I_{sat}$')
    ax[0].set_ylabel(r'Line-centre absorption $α$')
    h, _ = ax[1].get_legend_handles_labels()
    ax[1].legend(ncol=4, handles=[lot, lct, *h], frameon=False, loc='lower center', bbox_to_anchor=(0.5, 1.0))
    plt.margins(x=0)
    plt.tight_layout()
    fig.subplots_adjust(hspace=0, wspace=0)
    plt.savefig(f'fig4.png', dpi=300)
    plt.show()
