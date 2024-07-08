# PyMCSP. Python and Machine Learning for Crystal Structure Prediction and Diffraction Study

One of the most important and critical problems in materials science is determining the crystal structure or structures of a given material. The positions that ions occupy in the unit cell of a crystal, along with their placement in the periodic table of elements, fully determine the physical, chemical, and functional properties of the material. These properties range from whether the material is an insulator, conductor, or semiconductor to more exotic phenomena such as superconductivity. Thus, it is essential to know the crystal structures of different materials. Scattering experiments, such as X-ray diffraction, help experimental scientists determine these structures. However, sometimes we are interested in performing simulations of materials that have not yet been synthesized and about which we know nothing of their structural properties.

Crystal structure prediction methods [[1]](#1) are methodologies that search for local (and hopefully global) minima of energy to find the stable and metastable phases of crystal materials. These methods usually perform well and, in some cases, are able to discover novel phases. However, they require the use of supercomputer facilities and a significant amount of computer time, as they typically involve many density functional theory (DFT) simulations.

In this context, we present **PyMCSP**, a **P**ython implementation and **M**achine-Learning Interatomic Potentials for **C**rystal **S**tructure **P**rediction and Diffraction Stuyd. This method is intended to find novel phases for well-known or new materials within reasonable time frames, ranging from a few minutes to a few hours on a standard user's computer.

## How it works

Given a chemical compound and its stoichiometry, the program constructs random crystal structures with phase groups compatible with the given stoichiometry. To achieve this, the [PyXtal](https://github.com/qzhu2017/PyXtal) [[2]](#2) Python library is used.

For the generated structures, ionic relaxation is performed using [M3GNet](https://github.com/materialsvirtuallab/m3gnet) [[3]](#3). M3GNet implements machine-learning interatomic potentials, allowing relaxations to be performed in a few seconds on a regular computer (instead of minutes or hours on a supercomputer). After this, it is strongly recommended to perform DFT relaxations on the most interesting resulting phases.

The relaxed structures are interesting enough to conclude the process, but it is also possible to distort the structures (with Gaussian noise) and relax them again to explore unexplored low-energy phases. This is implemented using a simple Metropolis algorithm: the ions of the lowest energy phases after the initial relaxation are distorted, and the new phases are relaxed. After these relaxations, their energies are computed. If the new phase has a lower energy than the original, it is accepted; otherwise, it is rejected. Repeating this process enough times (in what we call the Generations Loop) should allow us to explore the unreached minima in the energy hypersurface.

The program also allows for the search of crystal phases under fixed pressure. The desired pressure can be input, and after relaxing the found structures, the program will determine the compressed structures and energies for the given pressure. Currently, the program only accounts for compressions or expansions of the unit cell. Given a pressure, the volume of the distorted structure can be found by minimizing the enthalpy, $H$, which can be computed as:

$H(p,V)=E(V)+pV,$

where $V$ has to be selected in a way that minimizes $H$ for a fixed pressure:

$\left. \frac{\partial H(p,V)}{\partial V} \right|_{p} = 0.$

After this, the pressurized structures are ranked based on enthalpy (instead of internal energy, as in computations without pressure).

Since M3GNet can be retrained with DFT results, we can use DFT materials simulation results to retrain M3GNet to look for novel phases of similar materials.

Using PyMCSP, we can predict the phase group of a given experimental diffractogram. To do this, a CSV file with the intensities as a function of $2\theta$ should be provided. The diffraction script will compare the experimental curves $I_{exp}$ with the theoretical curves $I_{theo}$ by computing a loss factor:

$Loss=a\left|n_{exp} - n_{theo} \right|^{2} + \int_{2\theta_{min}}^{2\theta_{max}}\left| I_{exp} -I_{theo} \right|^{2}d2\theta,$

here, the first term penalizes large differences in peak numbers, while the second term penalizes peaks in different $2\theta$ positions. After the program executes, the results are displayed in a file. The program also slightly varies the lattice parameters and selects the volume that minimizes the loss factor, thereby matching the experimental lattice parameters.

## Requirements

The code has been tested on the GNU/Linux operating system (Ubuntu 22.04.3 LTS) with Python 3.10.12. The required packages to execute PyMCSP are: 
- numpy
- scipy
- pymatgen
- pyxtal
- matgl 1.1.2

These packages can be downloaded manually, or you can execute the following command to install them automatically:
```bash
$ pip install -r requirements.txt
```

## Installation

To download the repository, use:

```bash
$ git clone https://github.com/polbeni/PyMCSP
```

Make sure to install all the necessary packages.

## How to use it

For any calculation, execute the main script `src/PyMCSP.py` with:

```bash
$ cd src
$ python3 PyMCSP.py
```

#### Crystal Structure Prediction

To perform crystal structure prediction, select the option `1` and verify that all the inputs are correct. The calculations will start, and a directory named `structures_files` will store the generated and relaxed structures. If you want to perform calculations with different inputs (for example, with different materials), you should save the `structures_files` directory in another location, as each time PyMCSP is executed, it overwrites this directory.

Pressure computations are performed by selecting the option `2`. The results will be saved in `structures_files/pressure_structures`. If you are interested in performing pressure calculations at different pressures, please save the `structures_files/pressure_structures` directory in a different location; otherwise, the previous results will be overwritten.

To apply the generations loop implementation, select the option `3`. The resulting structures will be stored in `structures_files/generations`. The results of generation number $i$ will be saved in `structures_files/generations/generation-i`, while the final structures will be in `structures_files/generations/final_structures`.

To perform pressure calculations or generations loops, it is necessary to have previously performed a crystal structure prediction calculation with option `1`, or provide the path to some of these results in the input file.

#### Diffraction Study

To do diffraction studies, select the option `4`.A CSV file with the experimental data should be provided, with the $2\theta$ values in the first column and $I$ A crystal structure prediction calculation must have been performed previously, as the program uses the found structures to compare with the experimental diffraction data. From the experimental diffractogram, the software reconstructs a clean diffraction curve to compare with the theoretically computed curves. Using option `5`, it is possible to tune the generation of the clean diffraction curve. Manually introduce the number of prominence and width of the peaks to identify the peaks of interest.

## Authors

This code and repository are being developed by:
- Pol Benítez Colominas (pol.benitez@upc.edu)
- Claudio Cazorla Silva (claudio.cazorla@upc.edu)

## References

<a id="1">[1]</a> 
WOODLEY, Scott M.; CATLOW, Richard. Crystal structure prediction from first principles. <em>Nature materials</em>, 2008, 7.12: 937-946.

<a id="2">[2]</a> 
FREDERICKS, Scott, et al. PyXtal: A Python library for crystal structure generation and symmetry analysis. <em>Computer Physics Communications</em>, 2021, 261: 107810.

<a id="3">[3]</a> 
CHEN, Chi; ONG, Shyue Ping. A universal graph deep learning interatomic potential for the periodic table. <em>Nature Computational Science</em>, 2022, 2.11: 718-728.


