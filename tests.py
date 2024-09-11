import matplotlib.pyplot as plt
import numpy as np
import os
import scipy as sp
from datetime import datetime
from scipy import constants as c
import sys
import importlib
import pkgutil
import elecsus
def import_submodules(module):
    """Import all submodules of a module, recursively."""
    for loader, module_name, is_pkg in pkgutil.walk_packages(
            module.__path__, module.__name__ + '.'):
        importlib.import_module(module_name)
import_submodules(elecsus)
import LindbladMasterEq as LME

package_name = 'Ahint'
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=plt.cm.tab20.colors)
results_folder = 'LME_plots'
os.makedirs(results_folder, exist_ok=True)
groundState = LME.state(5, 0, 1/2)	 # 5S1/2
excitedState_D1 = LME.state(5, 1, 1/2)	# 5P1/2
excitedState_D2 = LME.state(5, 1, 3/2)	# 5P3/2
E_LCP = elecsus.libs.BasisChanger.lrz_to_xyz(np.array([1,0,0]))
E_RCP = elecsus.libs.BasisChanger.lrz_to_xyz(np.array([0,1,0]))
E_LP = np.array([1,0,0])

###############################################################################
# Test all alkali atoms
###############################################################################
def sodium_default():
	p_dict_D1 = {'Elem': 'Na', 'Dline':'D1', 'lcell':2e-3, 'T': 20., 'laserPower': 1e-15, 'laserWaist': 5e-3}
	p_dict_D2 = {**p_dict_D1, 'Dline':'D2'}
	x = np.linspace(-4000, 4000, 2000)
	[y_elecsus_D1] = elecsus.elecsus_methods.calculate(x, E_in=E_LCP, p_dict=p_dict_D1, outputs=['S0'])
	[y_elecsus_D2] = elecsus.elecsus_methods.calculate(x, E_in=E_LCP, p_dict=p_dict_D2, outputs=['S0'])
	y_bwf_D1 = LME.get_spectra(x, E_in=E_LP, p_dict=p_dict_D1)
	y_bwf_D2 = LME.get_spectra(x, E_in=E_LP, p_dict=p_dict_D2)
	plt.figure()
	plt.plot(x, y_bwf_D1, c='C1', label='LME D1')
	plt.plot(x, y_elecsus_D1, '--', c='C0', label='ElecSus D1')
	plt.plot(x, y_bwf_D2, c='C3', label='LME D2')
	plt.plot(x, y_elecsus_D2, '--', c='C2', label='ElecSus D2')
	plt.xlabel('Detuning [MHz]')
	plt.ylabel('Transmission')
	plt.legend(frameon=False)
	plt.savefig(f'{results_folder}/sodium.png', dpi=200)

def potassium_default():
	p_dict_D1 = {'Elem': 'K', 'Dline':'D1', 'lcell':2e-3, 'T': 20., 'laserPower': 1e-15, 'laserWaist': 5e-3}
	p_dict_D2 = {**p_dict_D1, 'Dline':'D2'}
	x = np.linspace(-2000, 2000, 1000)
	[y_elecsus_D1] = elecsus.elecsus_methods.calculate(x, E_in=E_LCP, p_dict=p_dict_D1, outputs=['S0'])
	[y_elecsus_D2] = elecsus.elecsus_methods.calculate(x, E_in=E_LCP, p_dict=p_dict_D2, outputs=['S0'])
	y_bwf_D1 = LME.get_spectra(x, E_in=E_LP, p_dict=p_dict_D1)
	y_bwf_D2 = LME.get_spectra(x, E_in=E_LP, p_dict=p_dict_D2)
	plt.figure()
	plt.plot(x, y_bwf_D1, c='C1', label='LME D1')
	plt.plot(x, y_elecsus_D1, '--', c='C0', label='ElecSus D1')
	plt.plot(x, y_bwf_D2, c='C3', label='LME D2')
	plt.plot(x, y_elecsus_D2, '--', c='C2', label='ElecSus D2')
	plt.xlabel('Detuning [MHz]')
	plt.ylabel('Transmission')
	plt.legend(frameon=False)
	plt.savefig(f'{results_folder}/potassium.png', dpi=200)

def rubidium_default():
	p_dict_D1 = {'Elem': 'Rb', 'Dline':'D1', 'lcell':2e-3, 'T': 20., 'laserPower': 1e-15, 'laserWaist': 5e-3}
	p_dict_D2 = {**p_dict_D1, 'Dline':'D2'}
	x = np.linspace(-4000, 6000, 2000)
	[y_elecsus_D1] = elecsus.elecsus_methods.calculate(x, E_in=E_LCP, p_dict=p_dict_D1, outputs=['S0'])
	[y_elecsus_D2] = elecsus.elecsus_methods.calculate(x, E_in=E_LCP, p_dict=p_dict_D2, outputs=['S0'])
	y_bwf_D1 = LME.get_spectra(x, E_in=E_LP, p_dict=p_dict_D1)
	y_bwf_D2 = LME.get_spectra(x, E_in=E_LP, p_dict=p_dict_D2)
	plt.figure()
	plt.plot(x, y_bwf_D1, c='C1', label='LME D1')
	plt.plot(x, y_elecsus_D1, '--', c='C0', label='ElecSus D1')
	plt.plot(x, y_bwf_D2, c='C3', label='LME D2')
	plt.plot(x, y_elecsus_D2, '--', c='C2', label='ElecSus D2')
	plt.xlabel('Detuning [MHz]')
	plt.ylabel('Transmission')
	plt.legend(frameon=False)
	plt.savefig(f'{results_folder}/rubidium.png', dpi=200)

def caesium_default():
	p_dict_D1 = {'Elem': 'Cs', 'Dline':'D1', 'lcell':2e-3, 'T': 20., 'laserPower': 1e-15, 'laserWaist': 5e-3}
	p_dict_D2 = {**p_dict_D1, 'Dline':'D2'}
	x = np.linspace(-4000, 6000, 1000)
	[y_elecsus_D1] = elecsus.elecsus_methods.calculate(x, E_in=E_LCP, p_dict=p_dict_D1, outputs=['S0'])
	[y_elecsus_D2] = elecsus.elecsus_methods.calculate(x, E_in=E_LCP, p_dict=p_dict_D2, outputs=['S0'])
	y_bwf_D1 = LME.get_spectra(x, E_in=E_LP, p_dict=p_dict_D1)
	y_bwf_D2 = LME.get_spectra(x, E_in=E_LP, p_dict=p_dict_D2)
	plt.figure()
	plt.plot(x, y_bwf_D1, c='C1', label='LME D1')
	plt.plot(x, y_elecsus_D1, '--', c='C0', label='ElecSus D1')
	plt.plot(x, y_bwf_D2, c='C3', label='LME D2')
	plt.plot(x, y_elecsus_D2, '--', c='C2', label='ElecSus D2')
	plt.xlabel('Detuning [MHz]')
	plt.ylabel('Transmission')
	plt.legend(frameon=False)
	plt.savefig(f'{results_folder}/caesium.png', dpi=200)


###############################################################################
# Rb85 D1
###############################################################################
########## LCP ##########
def Rb85_D1_LCP_B_0G():
	p_dict = {'Elem':'Rb','Dline':'D1', 'lcell':2e-3, 'T': 20.,
	   'Bfield': 0, 'rb85frac': 100, 'Constrain': False, 'DoppTemp': -273.14999,
	   'laserPower': 1e-15, 'laserWaist': 5e-3}
	x = np.linspace(1400, 2100, 1000)
	[y_elecsus] = elecsus.elecsus_methods.calculate(x, E_in=E_LCP, p_dict=p_dict, outputs=['S0'])
	y_bwf = LME.get_spectra(x, E_in=E_LCP, p_dict=p_dict)
	plt.figure()
	plt.plot(x, y_bwf, c='C1', label=package_name)
	plt.plot(x, y_elecsus, '--', c='C4', label='ElecSus')
	plt.xlabel('Detuning [MHz]')
	plt.ylabel('Transmission')
	plt.legend(frameon=False)
	plt.savefig(f'{results_folder}/Rb85_D1_LCP_B_0G.png', dpi=200)

def Rb85_D1_LCP_B_100G():
	p_dict = {'Elem':'Rb','Dline':'D1', 'lcell':2e-3, 'T': 20.,
	   'Bfield': 100, 'rb85frac': 100, 'Constrain': False, 'DoppTemp': -273.14999, 'symbolic_transit': True,
	   'laserPower': 1e-15, 'laserWaist': 2e-3}
	x = np.linspace(1400, 2100, 1000)
	[y_elecsus] = elecsus.elecsus_methods.calculate(x, E_in=E_LCP, p_dict=p_dict, outputs=['S0'])
	y_bwf = LME.get_spectra(x, E_in=E_LCP, p_dict=p_dict)
	fig = plt.figure(tight_layout=True)
	plt.plot(x, y_bwf, c='C1', label=package_name)
	plt.plot(x, y_elecsus, '--', c='C4', label='ElecSus')
	plt.xlabel('Detuning [MHz]')
	plt.ylabel('Transmission')
	plt.legend(frameon=False)
	plt.savefig(f'{results_folder}/Rb85_D1_LCP_B_100G.png', dpi=200)

def Rb85_D1_LCP_B_1000G():
	p_dict = {'Elem':'Rb','Dline':'D1', 'lcell':2e-3, 'T': 20.,
	   'Bfield': 1000, 'rb85frac': 100, 'Constrain': False, 'DoppTemp': -273.14999,
		'laserPower': 1e-15, 'laserWaist': 5e-3}
	x = np.linspace(1400, 2100, 1000)
	[y_elecsus] = elecsus.elecsus_methods.calculate(x, E_in=E_LCP, p_dict=p_dict, outputs=['S0'])
	y_bwf = LME.get_spectra(x, E_in=E_LCP, p_dict=p_dict)
	plt.figure()
	plt.plot(x, y_bwf, c='C1', label=package_name)
	plt.plot(x, y_elecsus, '--', c='C4', label='ElecSus')
	plt.xlabel('Detuning [MHz]')
	plt.ylabel('Transmission')
	plt.legend(frameon=False)
	plt.savefig(f'{results_folder}/Rb85_D1_LCP_B_1000G.png', dpi=200)

########## RCP ##########
def Rb85_D1_RCP_B_0G():
	p_dict = {'Elem':'Rb','Dline':'D1', 'lcell':2e-3, 'T': 20.,
	   'Bfield': 0, 'rb85frac': 100, 'Constrain': False, 'DoppTemp': -273.14999,
	   'laserPower': 1e-15, 'laserWaist': 5e-3}
	x = np.linspace(1400, 2100, 1000)
	[y_elecsus] = elecsus.elecsus_methods.calculate(x, E_in=E_RCP, p_dict=p_dict, outputs=['S0'])
	y_bwf = LME.get_spectra(x, E_in=E_RCP, p_dict=p_dict)
	plt.figure()
	plt.plot(x, y_bwf, c='C1', label=package_name)
	plt.plot(x, y_elecsus, '--', c='C4', label='ElecSus')
	plt.xlabel('Detuning [MHz]')
	plt.ylabel('Transmission')
	plt.legend(frameon=False)
	plt.savefig(f'{results_folder}/Rb85_D1_RCP_B_0G.png', dpi=200)

def Rb85_D1_RCP_B_100G():
	p_dict = {'Elem':'Rb','Dline':'D1', 'lcell':2e-3, 'T': 20.,
	   'Bfield': 100, 'rb85frac': 100, 'Constrain': False, 'DoppTemp': -273.14999,
	   'laserPower': 1e-15, 'laserWaist': 5e-3}
	x = np.linspace(1400, 2100, 1000)
	[y_elecsus] = elecsus.elecsus_methods.calculate(x, E_in=E_RCP, p_dict=p_dict, outputs=['S0'])
	y_bwf = LME.get_spectra(x, E_in=E_RCP, p_dict=p_dict)
	plt.figure()
	plt.plot(x, y_bwf, c='C1', label=package_name)
	plt.plot(x, y_elecsus, '--', c='C4', label='ElecSus')
	plt.xlabel('Detuning [MHz]')
	plt.ylabel('Transmission')
	plt.legend(frameon=False)
	plt.savefig(f'{results_folder}/Rb85_D1_RCP_B_100G.png', dpi=200)

def Rb85_D1_RCP_B_1000G():
	p_dict = {'Elem':'Rb','Dline':'D1', 'lcell':2e-3, 'T': 20.,
	   'Bfield': 1000, 'rb85frac': 100, 'Constrain': False, 'DoppTemp': -273.14999,
	   'laserPower': 1e-15, 'laserWaist': 5e-3}
	x = np.linspace(1400, 2100, 1000)
	[y_elecsus] = elecsus.elecsus_methods.calculate(x, E_in=E_RCP, p_dict=p_dict, outputs=['S0'])
	y_bwf = LME.get_spectra(x, E_in=E_RCP, p_dict=p_dict)
	plt.figure()
	plt.plot(x, y_bwf, c='C1', label=package_name)
	plt.plot(x, y_elecsus, '--', c='C4', label='ElecSus')
	plt.xlabel('Detuning [MHz]')
	plt.ylabel('Transmission')
	plt.legend(frameon=False)
	plt.savefig(f'{results_folder}/Rb85_D1_RCP_B_1000G.png', dpi=200)
###############################################################################
# Rb85 D2
###############################################################################
########## LCP ##########
def Rb85_D2_LCP_B_0G():
	p_dict = {'Elem':'Rb','Dline':'D2', 'lcell':2e-3, 'T': 20.,
	   'Bfield': 0, 'rb85frac': 100, 'Constrain': False, 'DoppTemp': -273.14999,
	   'laserPower': 1e-15, 'laserWaist': 5e-3}
	x = np.linspace(1400, 2100, 1000)
	[y_elecsus] = elecsus.elecsus_methods.calculate(x, E_in=E_LCP, p_dict=p_dict, outputs=['S0'])
	y_bwf = LME.get_spectra(x, E_in=E_LCP, p_dict=p_dict)
	plt.figure()
	plt.plot(x, y_bwf, c='C1', label=package_name)
	plt.plot(x, y_elecsus, '--', c='C4', label='ElecSus')
	plt.xlabel('Detuning [MHz]')
	plt.ylabel('Transmission')
	plt.legend(frameon=False)
	plt.savefig(f'{results_folder}/Rb85_D2_LCP_B_0G.png', dpi=200)

def Rb85_D2_LCP_B_100G():
	p_dict = {'Elem':'Rb','Dline':'D2', 'lcell':2e-3, 'T': 20.,
	   'Bfield': 100, 'rb85frac': 100, 'Constrain': False, 'DoppTemp': -273.14999,
	   'laserPower': 1e-15, 'laserWaist': 5e-3}
	x = np.linspace(1400, 2100, 1000)
	[y_elecsus] = elecsus.elecsus_methods.calculate(x, E_in=E_LCP, p_dict=p_dict, outputs=['S0'])
	y_bwf = LME.get_spectra(x, E_in=E_LCP, p_dict=p_dict)
	plt.figure()
	plt.plot(x, y_bwf, c='C1', label=package_name)
	plt.plot(x, y_elecsus, '--', c='C4', label='ElecSus')
	plt.xlabel('Detuning [MHz]')
	plt.ylabel('Transmission')
	plt.legend(frameon=False)
	plt.savefig(f'{results_folder}/Rb85_D2_LCP_B_100G.png', dpi=200)

########## RCP ##########
def Rb85_D2_RCP_B_100G():
	p_dict = {'Elem':'Rb','Dline':'D2', 'lcell':2e-3, 'T': 20.,
	   'Bfield': 100, 'rb85frac': 100, 'Constrain': False, 'DoppTemp': -273.14999,
	   'laserPower': 1e-15, 'laserWaist': 5e-3}
	x = np.linspace(1400, 2100, 1000)
	[y_elecsus] = elecsus.elecsus_methods.calculate(x, E_in=E_RCP, p_dict=p_dict, outputs=['S0'])
	y_bwf = LME.get_spectra(x, E_in=E_RCP, p_dict=p_dict)
	plt.figure()
	plt.plot(x, y_bwf, c='C1', label=package_name)
	plt.plot(x, y_elecsus, '--', c='C4', label='ElecSus')
	plt.xlabel('Detuning [MHz]')
	plt.ylabel('Transmission')
	plt.legend(frameon=False)
	plt.savefig(f'{results_folder}/Rb85_D2_RCP_B_100G.png', dpi=200)

###############################################################################
# Rb87 D1
###############################################################################
########## LCP ##########
def Rb87_D1_LCP_B_0G():
	p_dict = {'Elem':'Rb','Dline':'D1', 'lcell':2e-3, 'T': 20.,
	   'Bfield': 0, 'rb85frac': 0, 'Constrain': False, 'DoppTemp': -273.14999,
	   'laserPower': 1e-15, 'laserWaist': 5e-3}
	x = np.linspace(3500, 5000, 1000)
	[y_elecsus] = elecsus.elecsus_methods.calculate(x, E_in=E_LCP, p_dict=p_dict, outputs=['S0'])
	y_bwf = LME.get_spectra(x, E_in=E_LCP, p_dict=p_dict)
	plt.figure()
	plt.plot(x, y_bwf, c='C1', label=package_name)
	plt.plot(x, y_elecsus, '--', c='C4', label='ElecSus')
	plt.xlabel('Detuning [MHz]')
	plt.ylabel('Transmission')
	plt.legend(frameon=False)
	plt.savefig(f'{results_folder}/Rb87_D1_LCP_B_0G.png', dpi=200)

def Rb87_D1_LCP_B_100G():
	p_dict = {'Elem':'Rb', 'Dline':'D1', 'lcell':2e-3, 'T': 20.,
	   'Bfield': 100, 'rb85frac': 0, 'Constrain': False, 'DoppTemp': -273.14999,
	   'laserPower': 1e-15, 'laserWaist': 5e-3}
	x = np.linspace(4500, 4750, 1000)
	[y_elecsus] = elecsus.elecsus_methods.calculate(x, E_in=E_LCP, p_dict=p_dict, outputs=['S0'])
	y_bwf = LME.get_spectra(x, E_in=E_LCP, p_dict=p_dict)
	plt.figure()
	plt.plot(x, y_bwf, c='C1', label=package_name)
	plt.plot(x, y_elecsus, '--', c='C4', label='ElecSus')
	plt.xlabel('Detuning [MHz]')
	plt.ylabel('Transmission')
	plt.legend(frameon=False)
	plt.savefig(f'{results_folder}/Rb87_D1_LCP_B_100G.png', dpi=200)

def Rb87_D1_LCP_B_1000G():
	p_dict = {'Elem':'Rb','Dline':'D1', 'lcell':2e-3, 'T': 20.,
	   'Bfield': 1000, 'rb85frac': 0, 'Constrain': False, 'DoppTemp': -273.14999,
	   'laserPower': 1e-15, 'laserWaist': 5e-3}
	x = np.linspace(4250, 6000, 1000)
	[y_elecsus] = elecsus.elecsus_methods.calculate(x, E_in=E_LCP, p_dict=p_dict, outputs=['S0'])
	y_bwf = LME.get_spectra(x, E_in=E_LCP, p_dict=p_dict)
	plt.figure()
	plt.plot(x, y_bwf, c='C1', label=package_name)
	plt.plot(x, y_elecsus, '--', c='C4', label='ElecSus')
	plt.xlabel('Detuning [MHz]')
	plt.ylabel('Transmission')
	plt.legend(frameon=False)
	plt.savefig(f'{results_folder}/Rb87_D1_LCP_B_1000G.png', dpi=200)

def Rb87_D1_LCP_B_6000G():
	p_dict = {'Elem':'Rb','Dline':'D1', 'lcell':2e-3, 'T': 20.,
	   'Bfield': 6000, 'rb85frac': 0, 'Constrain': False, 'DoppTemp': -273.14999,
	   'laserPower': 1e-15, 'laserWaist': 5e-3}
	x = np.linspace(5500, 9000, 1000)
	[y_elecsus] = elecsus.elecsus_methods.calculate(x, E_in=E_LCP, p_dict=p_dict, outputs=['S0'])
	y_bwf = LME.get_spectra(x, E_in=E_LCP, p_dict=p_dict)
	plt.figure()
	plt.plot(x, y_bwf, c='C1', label=package_name)
	plt.plot(x, y_elecsus, '--', c='C4', label='ElecSus')
	plt.xlabel('Detuning [MHz]')
	plt.ylabel('Transmission')
	plt.legend(frameon=False)
	plt.savefig(f'{results_folder}/Rb87_D1_LCP_B_6000G.png', dpi=200)

########## RCP ##########
def Rb87_D1_RCP_B_0G():
	p_dict = {'Elem':'Rb','Dline':'D1', 'lcell':2e-3, 'T': 20.,
	   'Bfield': 0, 'rb85frac': 0, 'Constrain': False, 'DoppTemp': -273.14999,
	   'laserPower': 1e-15, 'laserWaist': 5e-3}
	x = np.linspace(3500, 5000, 1000)
	[y_elecsus] = elecsus.elecsus_methods.calculate(x, E_in=E_RCP, p_dict=p_dict, outputs=['S0'])
	y_bwf = LME.get_spectra(x, E_in=E_RCP, p_dict=p_dict)
	plt.figure()
	plt.plot(x, y_bwf, c='C1', label=package_name)
	plt.plot(x, y_elecsus, '--', c='C4', label='ElecSus')
	plt.xlabel('Detuning [MHz]')
	plt.ylabel('Transmission')
	plt.legend(frameon=False)
	plt.savefig(f'{results_folder}/Rb87_D1_RCP_B_0G.png', dpi=200)

def Rb87_D1_RCP_B_100G():
	p_dict = {'Elem':'Rb','Dline':'D1', 'lcell':2e-3, 'T': 20.,
	   'Bfield': 100, 'rb85frac': 0, 'Constrain': False, 'DoppTemp': -273.14999,
	   'laserPower': 1e-15, 'laserWaist': 5e-3}
	x = np.linspace(3500, 5000, 2000)
	[y_elecsus] = elecsus.elecsus_methods.calculate(x, E_in=E_RCP, p_dict=p_dict, outputs=['S0'])
	y_bwf = LME.get_spectra(x, E_in=E_RCP, p_dict=p_dict)
	plt.figure()
	plt.plot(x, y_bwf, c='C1', label=package_name)
	plt.plot(x, y_elecsus, '--', c='C4', label='ElecSus')
	plt.xlabel('Detuning [MHz]')
	plt.ylabel('Transmission')
	plt.legend(frameon=False)
	plt.savefig(f'{results_folder}/Rb87_D1_RCP_B_100G.png', dpi=200)

def Rb87_D1_RCP_B_1000G():
	p_dict = {'Elem':'Rb','Dline':'D1', 'lcell':2e-3, 'T': 20.,
	   'Bfield': 1000, 'rb85frac': 0, 'Constrain': False, 'DoppTemp': -273.14999,
	   'laserPower': 1e-15, 'laserWaist': 5e-3}
	x = np.linspace(3500, 6000, 1000)
	[y_elecsus] = elecsus.elecsus_methods.calculate(x, E_in=E_RCP, p_dict=p_dict, outputs=['S0'])
	y_bwf = LME.get_spectra(x, E_in=E_RCP, p_dict=p_dict)
	plt.figure()
	plt.plot(x, y_bwf, c='C1', label=package_name)
	plt.plot(x, y_elecsus, '--', c='C4', label='ElecSus')
	plt.xlabel('Detuning [MHz]')
	plt.ylabel('Transmission')
	plt.legend(frameon=False)
	plt.savefig(f'{results_folder}/Rb87_D1_RCP_B_1000G.png', dpi=200)

def Rb87_D1_RCP_B_6000G():
	p_dict = {'Elem':'Rb','Dline':'D1', 'lcell':2e-3, 'T': 20.,
	   'Bfield': 6000, 'rb85frac': 0, 'Constrain': False, 'DoppTemp': -273.14999,
	   'laserPower': 1e-15, 'laserWaist': 5e-3}
	x = np.linspace(8000, 15000, 10000)
	[y_elecsus] = elecsus.elecsus_methods.calculate(x, E_in=E_RCP, p_dict=p_dict, outputs=['S0'])
	y_bwf = LME.get_spectra(x, E_in=E_RCP, p_dict=p_dict)
	plt.figure()
	plt.plot(x, y_bwf, c='C1', label=package_name)
	plt.plot(x, y_elecsus, '--', c='C4', label='ElecSus')
	plt.xlabel('Detuning [MHz]')
	plt.ylabel('Transmission')
	plt.legend(frameon=False)
	plt.savefig(f'{results_folder}/Rb87_D1_RCP_B_6000G.png', dpi=200)

########## LP ##########
def Rb87_D1_LP_B_0G():
	p_dict = {'Elem':'Rb','Dline':'D1', 'lcell':2e-3, 'T': 20.,
	   'Bfield': 0, 'rb85frac': 0, 'Constrain': False, 'DoppTemp': -273.14999,
	   'laserPower': 1e-15, 'laserWaist': 5e-3}
	x = np.linspace(3500, 5000, 1000)
	[y_elecsus] = elecsus.elecsus_methods.calculate(x, E_in=E_LP, p_dict=p_dict, outputs=['S0'])
	y_bwf = LME.get_spectra(x, E_in=E_LP, p_dict=p_dict)
	plt.figure()
	plt.plot(x, y_bwf, c='C1', label=package_name)
	plt.plot(x, y_elecsus, '--', c='C4', label='ElecSus')
	plt.xlabel('Detuning [MHz]')
	plt.ylabel('Transmission')
	plt.legend(frameon=False)
	plt.savefig(f'{results_folder}/Rb87_D1_LP_B_0G.png', dpi=200)

def Rb87_D1_LP_B_100G():
	p_dict = {'Elem':'Rb','Dline':'D1', 'lcell':2e-3, 'T': 20.,
	   'Bfield': 100, 'rb85frac': 0, 'Constrain': False, 'DoppTemp': -273.14999,
	   'laserPower': 1e-15, 'laserWaist': 5e-3}
	x = np.linspace(3500, 5000, 1000)
	[y_elecsus] = elecsus.elecsus_methods.calculate(x, E_in=E_LP, p_dict=p_dict, outputs=['S0'])
	y_bwf = LME.get_spectra(x, E_in=E_LP, p_dict=p_dict)
	plt.figure()
	plt.plot(x, y_bwf, c='C1', label=package_name)
	plt.plot(x, y_elecsus, '--', c='C4', label='ElecSus')
	plt.xlabel('Detuning [MHz]')
	plt.ylabel('Transmission')
	plt.legend(frameon=False)
	plt.savefig(f'{results_folder}/Rb87_D1_LP_B_100G.png', dpi=200)

def Rb87_D2_LCP_B_0G():
	p_dict = {'Elem':'Rb','Dline':'D2', 'lcell':2e-3, 'T': 50.,
	   'Bfield': 0, 'rb85frac': 0, 'Constrain': False, 'DoppTemp': -273.14999,
	   'laserPower': 1e-15, 'laserWaist': 5e-3}
	x = np.linspace(3900, 4400, 1000)
	# x = np.linspace(-6000, 6000, 10000)
	[y_elecsus] = elecsus.elecsus_methods.calculate(x, E_in=E_LCP, p_dict=p_dict, outputs=['S0'])
	y_bwf = LME.get_spectra(x, E_in=E_LCP, p_dict=p_dict)
	plt.figure()
	plt.plot(x, y_bwf, c='C1', label=package_name)
	plt.plot(x, y_elecsus, '--', c='C4', label='ElecSus')
	plt.xlabel('Detuning [MHz]')
	plt.ylabel('Transmission')
	plt.legend(frameon=False)
	plt.savefig(f'{results_folder}/Rb87_D2_LCP_B_0G.png', dpi=200)

def Rb87_D2_LCP_B_100G():
	p_dict = {'Elem':'Rb','Dline':'D2', 'lcell':2e-3, 'T': 20.,
	   'Bfield': 100, 'rb85frac': 0, 'Constrain': False, 'DoppTemp': -273.14999, 'symbolic_transit': True,
	   'laserPower': 1e-15, 'laserWaist': 2e-3}
	x = np.linspace(3500, 4500, 1000)
	[y_elecsus] = elecsus.elecsus_methods.calculate(x, E_in=E_LCP, p_dict=p_dict, outputs=['S0'])
	y_bwf = LME.get_spectra(x, E_in=E_LCP, p_dict=p_dict)
	plt.figure()
	plt.plot(x, y_bwf, c='C1', label=package_name)
	plt.plot(x, y_elecsus, '--', c='C4', label='ElecSus')
	plt.xlabel('Detuning [MHz]')
	plt.ylabel('Transmission')
	plt.legend(frameon=False)
	plt.savefig(f'{results_folder}/Rb87_D2_LCP_B_100G.png', dpi=200)


def Rb87_D1_LCP_B_100G_power_scan():
	p_dict = {'Elem':'Rb','Dline':'D1', 'lcell':2e-3, 'T': 10.,
	   'Bfield': 100, 'rb85frac': 0, 'Constrain': False, 'DoppTemp': -273.1499,
	   'laserPower': 1e-15, 'laserWaist': 5e-3}
	x = np.linspace(4500, 4550, 1000)
	[y_elecsus] = elecsus.elecsus_methods.calculate(x, E_in=E_LCP, p_dict=p_dict, outputs=['S0'])
	plt.figure()
	for p in [1e-10, 1e-8, 1e-7]:
		p_dict['laserPower'] = p
		y_bwf = LME.get_spectra(x, E_in=E_LCP, p_dict=p_dict)
		plt.plot(x, y_bwf, c='C1', label=package_name)
	plt.plot(x, y_elecsus, '--', c='C4', label='ElecSus')
	plt.xlabel('Detuning [MHz]')
	plt.ylabel('Transmission')
	plt.legend(frameon=False)
	plt.savefig(f'{results_folder}/Rb87_D1_LCP_B_100G_power_scan.png', dpi=200)

def Rb87_D1_LCP_B_100G_high_T():
	p_dict = {'Elem':'Rb','Dline':'D1', 'lcell':2e-3, 'T': 20., 'Bfield': 100, 'rb85frac': 0,
		   'laserPower': 1e-15, 'laserWaist': 5e-3}
	x = np.linspace(3000, 6000, 1000)
	[y_elecsus] = elecsus.elecsus_methods.calculate(x, E_in=E_LCP, p_dict=p_dict, outputs=['S0'])
	y_bwf = LME.get_spectra(x, E_in=E_LCP, p_dict=p_dict)
	plt.figure()
	plt.plot(x, y_bwf, c='C1', label=package_name)
	plt.plot(x, y_elecsus, '--', c='C4', label='ElecSus')
	plt.xlabel('Detuning [MHz]')
	plt.ylabel('Transmission')
	plt.legend(frameon=False)
	plt.savefig(f'{results_folder}/Rb87_D1_LCP_B_100G_high_T.png', dpi=200)

def Rb87_D2_RCP_B_6000G_high_T():
	from datetime import datetime
	p_dict = {'Elem':'Rb','Dline':'D2', 'lcell':2e-3, 'T': 64., 'Bfield': 6000, 'rb85frac': 0, 'GammaBuf': 0,
		   'laserPower': 1e-15, 'laserWaist': 2e-3, 'collisions': 'decay'}#, 'symbolic_transit': True}
	x = np.linspace(4000, 20000, 2000)
	[y_elecsus] = elecsus.elecsus_methods.calculate(x, E_in=E_RCP, p_dict=p_dict, outputs=['S0'])
	y_bwf = LME.get_spectra(x, E_in=E_RCP, p_dict=p_dict)
	fig = plt.figure(tight_layout=True)
	plt.plot(x, y_bwf, c='C1', label=package_name)
	plt.plot(x, y_elecsus, '--', c='C4', label='ElecSus')
	plt.xlabel('Detuning [MHz]')
	plt.ylabel('Transmission')
	plt.legend(frameon=False)
	plt.savefig(f'{results_folder}/Rb87_D2_RCP_B_6000G_high_T.png', dpi=200)

###############################################################################
# Collection of other tests
###############################################################################
def Rb87_D2_RCP_B_6000G_high_T_custom_transit():
	from datetime import datetime
	p_dict = {'Elem':'Rb','Dline':'D2', 'lcell':2e-3, 'T': 60., 'Bfield': 6000, 'rb85frac': 0, 'GammaBuf': 0}#, 'DoppTemp': -273.14999}
	p_dict_bwf = {**p_dict, 'laserPower': 1e-7, 'laserWaist': 2e-3, 'symbolic_transit': True}

	def maxwell_boltzmann(v):
		# https://en.wikipedia.org/wiki/Maxwell%E2%80%93Boltzmann_distribution
		return 4 * np.pi * np.sqrt(Rb87_D2.mass / (2 * c.pi * c.k * Rb87_D2.DoppT))**3 \
		* np.exp(-Rb87_D2.mass * v**2 / (2 * c.k * Rb87_D2.DoppT)) * v**2

	def rayleigh(v):
		# https://en.wikipedia.org/wiki/Rayleigh_distribution
		return 2 * np.pi * Rb87_D2.mass / (2 * c.pi * c.k * Rb87_D2.DoppT) \
		* np.exp(-Rb87_D2.mass * v**2 / (2 * c.k * Rb87_D2.DoppT)) * v

	# x = np.linspace(9500, 10900, 600)
	x = np.linspace(4000, 20000, 2000)
	t = datetime.now()
	Rb87_D2 = LME.atomicSystem('Rb87', [groundState, excitedState_D2], p_dict=p_dict_bwf)
	beam = LME.beam(w=x, P=1e-3, D=2e-3, profile='flat')
	doppler = True
	RbDen = Rb87_D2.getNumberDensity(Rb87_D2.T)
	k = Rb87_D2.f_resonance / c.c #/ 1e6

	v = np.linspace(0, 1000, 100000)
	v_dist_rayleigh = rayleigh(v)
	v_dist_maxwell_boltzmann = maxwell_boltzmann(v)
	v_avg_rayleigh = (v*v_dist_rayleigh).sum() / v_dist_rayleigh.sum()
	v_avg_maxwell_boltzmann = (v*v_dist_maxwell_boltzmann).sum() / v_dist_maxwell_boltzmann.sum()
	print(f'Average in 2D: {v_avg_rayleigh:.1f}')
	print(f'Average in 3D: {v_avg_maxwell_boltzmann:.1f}')
	plt.figure()
	plt.plot(v, v_dist_rayleigh, c='C1', label='Rayleigh distribution')
	plt.plot(v, v_dist_maxwell_boltzmann, '--', c='C4', label='MB distribution')
	plt.legend()
	plt.xlabel('Velocity [m/s]')
	plt.ylabel('P(v)')
	plt.savefig(f'{results_folder}/compare_velocity_distributions.png', dpi=200)

	Rb87_D2.update_transit(v_avg_rayleigh)
	y_avg_rayleigh = Rb87_D2.transmission([beam], doppler=doppler, z=2e-3, transit_type='single')
	y_rayleigh = Rb87_D2.transmission([beam], doppler=doppler, z=2e-3, transit_type='integral')
	print(f'Total time of call: {datetime.now()-t}')

	fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, tight_layout=True, gridspec_kw={'height_ratios': [3, 1]})
	ax1.plot(x, y_avg_rayleigh, c='C1', label='Rayleigh, avg')
	ax1.plot(x, y_rayleigh, '--', c='C4', label='Rayleigh, int')
	ax2.plot(x, (y_avg_rayleigh - y_rayleigh) * 100, c='C1')
	plt.xlabel('Detuning [MHz]')
	ax1.set_ylabel('Optical depth')
	ax2.set_ylabel('Residual [%]')
	ax1.legend(frameon=False)
	plt.savefig(f'{results_folder}/Compare_velocity_distributions_doppler_1e-07_power.png', dpi=200)

def Rb87_D2_RCP_B_6000G_high_T_custom_beam_shape():
	from datetime import datetime
	p_dict = {'Elem':'Rb','Dline':'D2', 'lcell':2e-3, 'T': 60., 'Bfield': 6000, 'rb85frac': 0, 'GammaBuf': 0}#, 'DoppTemp': -273.14999}
	p_dict_bwf = {**p_dict, 'laserPower': 1e-7, 'laserWaist': 2e-3, 'symbolic_transit': True}

	def clip(xy, z, D=2e-3):
		xv, yv = np.meshgrid(xy, xy)
		return np.where(np.sqrt(xv**2 + yv**2) < D/2, z, 0)

	w0_to_fwhm = lambda w0: np.sqrt(2*np.log(2)) * w0
	fwhm_to_e2 = lambda fwhm: 1.699 * fwhm
	e2_to_fwhm = lambda e2: e2 / 1.699
	fwhm_to_w0 = lambda fwhm: fwhm / np.sqrt(2*np.log(2))

	N = 501
	Nc = int((N - 1) / 2)
	P0 = 1e-4
	e2 = 4.0e-3
	w0 = fwhm_to_w0(e2_to_fwhm(e2))
	xy = np.linspace(-5e-3, 5e-3, N)
	xv, yv = np.meshgrid(xy, xy)
	gaussian = 2 * P0 / np.pi / w0**2 * np.exp(-2*(xv**2 + yv**2)/w0**2)
	gaussian_clipped = clip(xy, gaussian, 2e-3)
	print(f'Integral of Gaussian: {2 * np.pi * (gaussian[Nc, Nc:] * np.abs(xy[Nc:])).sum() * 10e-3}')

	flat = np.ones_like(gaussian) * gaussian.max()
	flat_clipped = clip(xy, flat, 2e-3)

	ptot_gaussian_clipped = (gaussian_clipped[Nc, Nc:] * np.abs(xy[Nc:])).sum()
	ptot_flat_clipped = (flat_clipped[Nc, Nc:] * np.abs(xy[Nc:])).sum()

	flat /= ptot_flat_clipped
	flat_clipped /= ptot_flat_clipped
	gaussian /= ptot_gaussian_clipped
	gaussian_clipped /= ptot_gaussian_clipped
	print(f'Integral of clipped flat top: {ptot_flat_clipped}')
	print(f'Integral of normalized clipped flat top: {(flat_clipped[Nc, Nc:] * np.abs(xy[Nc:])).sum()}')
	print((gaussian_clipped[Nc, Nc:] * np.abs(xy[Nc:])).sum())

	plt.figure()
	plt.pcolormesh(xy, xy, gaussian)
	plt.pcolormesh(xy, xy, gaussian_clipped)
	plt.figure()
	plt.pcolormesh(xy, xy, flat)
	plt.pcolormesh(xy, xy, flat_clipped)
	plt.figure()
	plt.plot(xy, gaussian[Nc], c='tab:green')
	plt.plot(xy, gaussian_clipped[Nc], c='tab:green', ls='--')
	plt.plot(xy, flat[Nc], c='tab:blue')
	plt.plot(xy, flat_clipped[Nc], c='tab:blue', ls='--')

	x = np.linspace(4000, 20000, 2000)
	# x = np.linspace(5000, 7000, 300)
	t = datetime.now()
	Rb87_D2 = LME.atomicSystem('Rb87', [groundState, excitedState_D2], p_dict=p_dict_bwf)
	RbDen = Rb87_D2.getNumberDensity(Rb87_D2.T)
	k = Rb87_D2.f_resonance / c.c
	print('d_xy')
	print((xy[1] - xy[0]))
	flat_tmp = 0
	gaussian_tmp = 0
	for i in range(Nc):
		print(flat_clipped[Nc, Nc+i])
		if flat_clipped[Nc, Nc+i] == 0:
			continue
		print(i)
		flat_beam = LME.beam(w=x, P=P0 * flat_clipped[Nc, Nc+i] / flat_clipped.max(), D=2e-3, profile='flat')
		gaussian_beam = LME.beam(w=x, P=P0 * gaussian_clipped[Nc, Nc+i] / flat_clipped.max(), D=2e-3, profile='flat')
		_, flat_chi = Rb87_D2.solve_w_doppler([flat_beam])
		_, gaussian_chi = Rb87_D2.solve_w_doppler([gaussian_beam])
		flat_tmp += flat_chi * flat_clipped[Nc, Nc+i] * xy[Nc+i]
		gaussian_tmp += gaussian_chi * gaussian_clipped[Nc, Nc+i] * xy[Nc+i]
	flat_integral = np.exp(4 * np.pi * k * np.sqrt(1.0 + flat_tmp * RbDen).imag * 2e-3)
	gaussian_integral = np.exp(4 * np.pi * k * np.sqrt(1.0 + gaussian_tmp * RbDen).imag * 2e-3)
	_, chi_avg = Rb87_D2.solve_w_doppler([LME.beam(w=x, P=P0, D=2e-3, profile='flat')])
	flat_avg = np.exp(4 * np.pi * k * np.sqrt(1.0 + chi_avg * RbDen).imag * 2e-3)


	print(f'Total time of call: {datetime.now()-t}')

	fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, tight_layout=True, gridspec_kw={'height_ratios': [3, 1]})
	ax1.plot(x, flat_avg, c='tab:red', label='Average')
	ax1.plot(x, flat_integral, '--', c='C4', label='Average flat')
	ax1.plot(x, gaussian_integral, '--', c='C2', label='Average gaussian')
	ax2.plot(x, (flat_integral - flat_avg) * 100, c='C4')
	ax2.plot(x, (gaussian_integral - flat_avg) * 100, c='C2')
	plt.xlabel('Detuning [MHz]')
	ax1.set_ylabel('Optical depth')
	ax2.set_ylabel('Residual [%]')
	ax1.legend(frameon=False)
	plt.savefig(f'{results_folder}/Compare_beam_shapes_1e-4.png', dpi=200)

def Rb87_D2_RCP_B_6000G_high_T_custom_transit_check_velocity_classes():
	from datetime import datetime
	p_dict = {'Elem':'Rb','Dline':'D2', 'lcell':2e-3, 'T': 50., 'Bfield': 6000, 'rb85frac': 0, 'GammaBuf': 0}#, 'DoppTemp': -273.14999}
	p_dict_bwf = {**p_dict, 'laserPower': 1e-15, 'laserWaist': 2e-3, 'collisions': 'decay', 'symbolic_transit': True}

	def rayleigh(v):
		# https://en.wikipedia.org/wiki/Rayleigh_distribution
		return 2 * np.pi * Rb87_D2.atom.mass / (2 * c.pi * c.k * Rb87_D2.DoppT) \
		* np.exp(-Rb87_D2.atom.mass * v**2 / (2 * c.k * Rb87_D2.DoppT)) * v

	x = np.linspace(10206.6, 10207, 2)
	t = datetime.now()
	Rb87_D2 = LME.atomicSystem('Rb87', [groundState, excitedState_D2], p_dict=p_dict_bwf)
	beam = LME.beam(w=x, P=1e-1, D=2e-3, profile='flat')
	doppler = True

	RbDen = Rb87_D2.atom.getNumberDensity(Rb87_D2.T)
	k = Rb87_D2.f_resonance / c.c #/ 1e6

	v = np.linspace(0, 1000, 100000)
	v_dist_rayleigh = rayleigh(v) * (v[1] - v[0])
	v_dist_rayleigh /= v_dist_rayleigh.sum()
	v_avg_rayleigh = (v*v_dist_rayleigh).sum() / v_dist_rayleigh.sum()
	print(f'Average in 2D: {v_avg_rayleigh:.1f}')

	Rb87_D2.update_transit(v_avg_rayleigh)
	y_avg_rayleigh = Rb87_D2.optical_depth([beam], doppler=doppler)

	N = np.arange(3, 41, dtype=int)
	y_rayleigh = np.zeros(N.size, dtype=np.complex128)

	for ni, n in enumerate(N):
		v = np.linspace(0, 900, n)
		dv = v[1] - v[0]
		tmp = 0
		for i, vi in enumerate(v):
			Rb87_D2.update_transit(vi)
			if not doppler:
				_, chi = Rb87_D2.solve([beam])
			else:
				_, chi = Rb87_D2.solve_w_doppler([beam])
			tmp += chi[0] * rayleigh(vi) * dv
		y_rayleigh[ni] = 4 * np.pi * k * np.sqrt(1.0 + tmp * RbDen).imag
		print(f'Total time of call: {datetime.now()-t}')

	plt.figure(tight_layout=True)
	plt.axhline(y_avg_rayleigh[0], c='k')
	plt.plot(N, y_rayleigh, c='C1', label='Rayleigh, int')
	plt.xlabel('Number of velocity samples')
	plt.ylabel('Line-centre optical depth')
	plt.savefig(f'{results_folder}/Convergence_number_velocity_classes.png', dpi=200)


if __name__ == '__main__':
	np.set_printoptions(linewidth=300)
	# sodium_default()
	# potassium_default()
	# rubidium_default()
	# caesium_default()

	# Rb85_D1_LCP_B_0G()
	# Rb85_D1_LCP_B_100G()
	# Rb85_D1_LCP_B_1000G()
	# Rb85_D1_RCP_B_0G()
	# Rb85_D1_RCP_B_100G()
	# Rb85_D1_RCP_B_1000G()
	# Rb85_D2_LCP_B_0G()
	# Rb85_D2_LCP_B_100G()
	# Rb85_D2_RCP_B_100G()

	# Rb87_D1_LCP_B_0G()
	Rb87_D1_LCP_B_100G()
	# Rb87_D1_LCP_B_1000G()
	# Rb87_D1_LCP_B_6000G()
	# Rb87_D1_RCP_B_0G()
	# Rb87_D1_RCP_B_100G()
	# Rb87_D1_RCP_B_1000G()
	# Rb87_D1_RCP_B_6000G()
	# Rb87_D1_LP_B_0G()
	# Rb87_D1_LP_B_100G()
	# Rb87_D2_LCP_B_0G()
	# Rb87_D2_LCP_B_100G()

	# Rb87_D1_LCP_B_100G_power_scan()
	# Rb87_D1_LCP_B_100G_high_T()
	# Rb87_D2_RCP_B_6000G_high_T()

	# Rb87_D2_RCP_B_6000G_high_T_custom_transit()
	# Rb87_D2_RCP_B_6000G_high_T_custom_beam_shape()
	# Rb87_D2_RCP_B_6000G_high_T_custom_transit_check_velocity_classes()

	plt.show()
