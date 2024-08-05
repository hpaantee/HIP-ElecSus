# Alkalis at HIgh iNTensities - Ahint

A program to calculate the complex-valued susceptibility of hot alkali vapour. Its particular advantage is the simulation of spectra at intensites above saturation. The program allows the calculation of the D1 and D2-line of sodium, potassium, rubidium and caesium, even in the presence of an magnetic field.

## Prerequisites
Must have installed a python 3 interpreter installed, with the following packages:
- ARC-alkali-rydberg-calculator 3.6
- ElecSus 3.0
- Matplotlib
- Numpy
- Scipy
- Symengine
- Sympy
- Importlib
- psutil
- Wxpython

## Usage
After installation of the required packages, examples are provided in the
`tests.py` file. To run them, un-comment the respective function in `main`.

The software provides two interfaces. The `get_spectra()` function provides a similar interface as the ElecSus software to calculate the estimated spectrum. Currently, most features are directly provided bythe `atomicSystem` class. This class allows to model the light-matter interaction for a single isotope.

A short example to calculate the D1 line of sodium would be
```python
# Define relevant atomic parameters
beamdiameter = 2e-3  # [m], to define the transit-time broadening
lcell = 2e-3  # [m]
p_dict = {'Elem': 'K', 'Dline':'D1', 'T': 20., 'Bfield': 100, 'lcell':lcell,
    'K40frac': 0, 'K41frac': 0, 'Constrain': False, 'DoppTemp': -273.1499,
    'laserWaist': beamdiameter}
E_in = np.array([1,0,0])  # linear polarization

detuning = np.linspace(-2000, 2500, 1500)  # [MHz]
# Define beam properties: detunings, power, diameter and beam profile (flat/gaussian)
laserPower = 1e-13  # [Watt]
laserbeam = beam(w=detuning, P=laserPower, D=beamdiameter, profile='gaussian')
isotope = atomicSystem('K39', E_in=E_in, p_dict=p_dict)
populations, susceptibility = isotope.solve([laserbeam])
transmission = isotope.transmission([laserbeam], doppler=False, z=lcell)
```

Note that when directly accessing the `atomicSystem` class, certain dictionary values are overwritten:
- `laserPower` and `laserWaist` are provided by the `beam` class; however, the laserWaist still needs to provided for the correct calculation of the transit-time broadening; it can also be changed retroactively by calling the `update_transit()` routine
- `atomicSystem` always corresponds to a single isotope with its natural abundance

## License
Apache like ElecSus?

## Change Log
v 1.0.0
- Initial release to the public