import glob
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
import numpy as np
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
import sys
import os
from elecsus import elecsus_methods as em

os.chdir('paper_figures')
E_in = np.array([np.sqrt(1/2),np.sqrt(1/2),0])

def fit(x, x0, lsx, T):
    p_dict = {'Elem': 'Rb', 'Dline': 'D2', 'T': T, 'lcell': 2e-3, 'shift':x0,
              'rb85frac': 0.70, 'Bfield':6000}
    f = em.spectra.get_spectra(x*lsx, E_in=E_in, p_dict=p_dict, outputs=['S0'])[0]
    # print(x.size, f.shape)
    # print(f)
    return f.real

def scale_x_axis(ref_spec, ref_baseline, cavity_spec):
    x, y_ref = np.loadtxt(ref_spec, delimiter=',', skiprows=5, unpack=True)
    _, y_cavity = np.loadtxt(cavity_spec, delimiter=',', skiprows=5, unpack=True)
    _, y_baseline = np.loadtxt(ref_baseline, delimiter=',', skiprows=5, unpack=True)
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    ax1.plot(x, y_cavity)
    ax2.plot(x, y_ref, c='C1')
    plt.title('Select range of interest')
    clip_borders = plt.ginput(n=2, timeout=0)
    clip_borders = sorted(clip_borders, key=lambda x: x[0])
    plt.close()
    clip_l = (np.abs(x - clip_borders[0][0])).argmin()
    clip_r = (np.abs(x - clip_borders[1][0])).argmin()

    x_clipped = x[clip_l:clip_r]
    # y_ref =  y_ref[clip_l:clip_r]
    y_cavity =  y_cavity[clip_l:clip_r]
    y_baseline =  y_baseline[clip_l:clip_r]
    # get frequency axis
    cavity_peaks, _ = find_peaks(-y_cavity, height=0.01, distance=300)
    f = np.arange(cavity_peaks.size) * FSR
    poly = np.poly1d(np.polyfit(x_clipped[cavity_peaks], f, deg=7))
    plt.figure()
    plt.plot(x, poly(x))
    plt.plot(x_clipped[cavity_peaks], f, 'x')
    plt.ginput(show_clicks=False)
    plt.close()

    plt.plot(x, y_ref)
    clip_borders = plt.ginput(n=2, timeout=0)
    clip_borders = sorted(clip_borders, key=lambda x: x[0])
    plt.close()
    clip_l = (np.abs(x - clip_borders[0][0])).argmin()
    clip_r = (np.abs(x - clip_borders[1][0])).argmin()
    xc = poly(x[clip_l:clip_r])
    y_ref =  y_ref[clip_l:clip_r]
    y_baseline =  y_baseline[clip_l:clip_r]
    plt.plot(xc, y_ref)
    s = plt.ginput(n=-1)
    plt.close()
    xc_bounds = xc[[(np.abs(xc - s[i][0])).argmin() for i in range(len(s))]]
    xc_fit = np.concatenate([xc[(xc>xc_bounds[i]) & (xc<xc_bounds[i+1])] for i in range(0, len(xc_bounds), 2)])
    y_cell_fit = np.concatenate([y_ref[(xc>xc_bounds[i]) & (xc<xc_bounds[i+1])] for i in range(0, len(xc_bounds), 2)])
    normalization_curve = np.poly1d(np.polyfit(xc_fit, y_cell_fit, deg=3))
    plt.figure()
    plt.plot(xc, y_ref)
    plt.plot(xc_fit, y_cell_fit)
    plt.plot(xc, normalization_curve(xc))
    plt.title('Select Rb87 left-most peak')
    s = plt.ginput(timeout=0)
    plt.close()
    y_ref /= normalization_curve(xc)
    shift = xc[(np.abs(xc - s[0][0])).argmin()]
    print(shift)
    xc -= shift - 6080

    p0 = [0, 1, 50]
    popt, _ = curve_fit(fit, xc, y_ref, p0=p0)

    plt.plot(xc, y_ref)
    plt.plot(xc, fit(xc, *popt))
    plt.plot(xc, fit(xc, *p0))
    plt.ginput(timeout=0, show_clicks=False)
    plt.close()
    print(popt[0])
    lsx = popt[1]
    print(f'lsx: {lsx}')
    print(popt)
    return lsx * xc, clip_borders, xc_bounds, clip_l, clip_r

def raw_processing(file_data, file_cavity, file_baseline):
    _, y_cell = np.loadtxt(file_data, delimiter=',', skiprows=5, unpack=True)
    _, y_cavity = np.loadtxt(file_cavity, delimiter=',', skiprows=5, unpack=True)
    _, y_baseline = np.loadtxt(file_baseline, delimiter=',', skiprows=5, unpack=True)
    y_cell =  y_cell[clip_l:clip_r]
    y_cavity = y_cavity[clip_l:clip_r]
    y_baseline =  y_baseline[clip_l:clip_r]
    base_curve = np.poly1d(np.polyfit(x, y_baseline, deg=3))
    y_cell -= base_curve(x)

    while True:
        class btnfoos:
            accept = False
            def redo(self, event):
                plt.close()
            def accept(self, event):
                self.accept = True
                plt.close()
        def update(val):
            deg = sdegree.val
            normalization_curve = np.poly1d(np.polyfit(xc_fit, y_cell_fit, deg=deg))
            l.set_ydata(normalization_curve(x))
            fig.canvas.draw_idle()

        fig, ax = plt.subplots()
        fig.subplots_adjust(bottom=0.2)
        ax.plot(x, y_cell)

        ax_degree = fig.add_axes([0.24, 0.05, 0.3, 0.075])
        ax_redo = fig.add_axes([0.59, 0.05, 0.1, 0.075])
        ax_accept = fig.add_axes([0.7, 0.05, 0.1, 0.075])
        sdegree = Slider(
            ax_degree, 'Poly deg', 1, 15,
            valinit=3, valstep = 1, initcolor=None
        )
        callback = btnfoos()
        sdegree.on_changed(update)
        bacc = Button(ax_accept, 'Accept')
        bacc.on_clicked(callback.accept)
        baredo = Button(ax_redo, 'Redo')
        baredo.on_clicked(callback.redo)

        ax.set_title('Click <= 1 for accepting spectrum')
        s = plt.ginput(n=-1)
        if len(s) <= 1:
            plt.close()
            break
        xc_bounds = x[[(np.abs(x - s[i][0])).argmin() for i in range(len(s))]]
        xc_fit = np.concatenate([x[(x>xc_bounds[i]) & (x<xc_bounds[i+1])] for i in range(0, len(xc_bounds), 2)])
        y_cell_fit = np.concatenate([y_cell[(x>xc_bounds[i]) & (x<xc_bounds[i+1])] for i in range(0, len(xc_bounds), 2)])
        ax.plot(xc_fit, y_cell_fit)
        normalization_curve = np.poly1d(np.polyfit(xc_fit, y_cell_fit, deg=3))
        l, = ax.plot(x, normalization_curve(x))
        plt.show()
        if callback.accept == True:
            deg = sdegree.val
            normalization_curve = np.poly1d(np.polyfit(xc_fit, y_cell_fit, deg=deg))
            y_cell /= normalization_curve(x)
    return x, y_cell

FSR = 360
dir = 'peak_measurements/2.0mm_0.6T_15V_bufgas'
prefix = 'poly7'

files = np.asarray(glob.glob(f'{dir}/C4--*W.csv'))
print(files)
powers = []
for s0 in files:
    s1 = s0.split('/')[-1]
    s2 = s1.split('--')[-1]
    s3 = s2.split('W')[0]
    val = float(s3[:-1])
    if s3[-1] == 'u':
        val *= 1e-6
    elif s3[-1] == 'm':
        val *= 1e-3
    else:
        raise ValueError('Wrong/missing unit!')
    powers.append(val)
powers = np.array(powers)
sort_idx = np.argsort(powers)
files = files[sort_idx]
powers = powers[sort_idx]

f_min = files[0]
print(f'Detected file {f_min} with minimum power')
f_min_baseline = f_min.split('.csv')[0] + '-baseline.csv'
f_min_cavity = f_min.replace("C4", "C2")

x, clip_borders, xc_bounds, clip_l, clip_r = scale_x_axis(f_min, f_min_baseline, f_min_cavity)
for s0 in files:
    s1 = s0.split('/')[-1]
    s2 = s1.split('--')[-1]
    s3 = s2.split('.csv')[0]
    if not os.path.isfile(f'{dir}/{prefix}/{s3}.npz'):
        print(f'Now processing file: {s3}')
        xt, yt = raw_processing(f'{dir}/C4--{s3}.csv', f'{dir}/C2--{s3}.csv', f_min_baseline)
        np.savez(f'{dir}/{prefix}/{s3}.npz', x=xt, y=yt)
