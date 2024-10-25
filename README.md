# High Intensity Probe ElecSus - HIP ElecSus

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
- tqdm
- Wxpython

## Usage
Poetry is used for the dependency management of packages. The required packages can either be installed manually or if poetry is available, by running the command  `poetry install`. This will create a virtual environment with all necessary packages, available via `poetry shell`.
After installation of the required packages, examples are provided in the
`tests.py` file. To run them, un-comment the respective function in `main`.
Alternatively, the examples are provided by jupyter notebooks, found in the `notebooks` folder.

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

## Performance
Solving the Lindblad master equation is computationally expensive and highly depends on the number of involved hyperfine states. This makes it rather slow. For orientation, below are listed typical times on a "normal" desktop for solving the linear equation system. All runs use a single isotope and 1000 values of detuning:
| Isotope | D line | Time [s] |
| :---: | :---: | ---: |
| Na   | D1 |   1 |
| Na   | D2 |   8 |
| K-39  | D1 |   1 |
| K-39  | D2 |   9 |
| Rb-85 | D1 |   8 |
| Rb-85 | D2 |  98 |
| Rb-87 | D1 |   1 |
| Rb-87 | D2 |   9 |
| Cs   | D1 |  53 |
| Cs   | D2 | 261 |
In terms of scaling, Rb87 D1 was scaled to increasing values of detuning:

| Detuning values | Time [s] |
| :---: | ---: |
| 1e3 |     1 |
| 1e4 |    17 |
| 1e5 |   181 |
| 1e6 |  1920 |

The values indicate a more or less linear scaling with the number of detuning values, as would be expected.

**Important note:** For Doppler broadening, the number of calculations is instead given by
$$ N = (\Delta f + 2000\,\text{MHz}) \cdot {0.5 \over \text{MHz}}, $$
where $\Delta f$ is the detuning range.

## Reproduction of Figure 3 and 4
The folder `paper_figures` contains all data necessary for reproducing the relevant figures of the paper. To run the respective scripts, call `python -m paper_figures.fig3`.
The subfolder `simulation` contains the line-centre absorptions simulated for Figure 4; `data` the experimentally measured data. `tmp_fig3` contains temporary files to speed up the figure generation.

## Shortcomings
Currently this software is not able to reproduce all the functionality of ElecSus. Especially, it does not allow
- arbitrary magnetic field orientations
- considering self-broadening at high temperatures
It also only allows computation of the transmission.
Note that this is not considered a todo list.

## License
All the files distributed with this program are provided subject to the Apache License, Version 2.0. A Copy of the license is provided.

## Change Log
v 1.0.0
- Initial release to the public