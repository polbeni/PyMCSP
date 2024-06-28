# PyMCSP. Python and Machine Learning for Crystal Structure Prediction and Diffraction Study

One of the most important and critical problems in materials science is to know the crystal structure or structures for a given material. The position that ions occupy in the unit cell of a crystal and in the periodic table of elements fully determines the physical, chemical and functional properties of the material, from if the material is an insulator, conductor or semiconductor to more exotic phenomena as superconductivity. Thus, it is important to know the crystal structures of the different materials. Scattering experiments (as X-ray diffraction) help experimental scientists to determine these structures, but sometimes we could be interested in performing simulations of materials that have not been synthesized yet and of which we know nothing of their structural properties. Crystal structure prediction methods [[1]](#1), are methodologies that look for local (and hopefully the global) minimums of energy, thus for the stable and metastable phases of crystal materials. These methods usually present a good performance and in some cases are able to find novel phases, but it is necessaty to use supercomputer facilities and a lot of computer time, since they usually perform plenty of DFT simulations.

In this context, we present **PyMCSP** a **P**ython and **M**achine Learning methods implementation for **C**rystal **S**tructure **P**rediction. This is intended to be a method to find novel phases for well known or new materials at reasonable frame times, few minutes to few hours in a standard user's computer.

The code in `src` is completely functional (old version) and it is possible to perform crystal structure prediction and diffraction studies. However, we are migrating to a more compact code and user friendly interface. The new code can be found in `src2`, although the diffraction part is still not available and it is not explained how to work with in this repository.

## How it works

From a given chemical compounds and stoichiometry, the program constructs random crystal structures with the phase groups that are compatible with the given stoichiometry. To do this [PyXtal](https://github.com/qzhu2017/PyXtal) [[2]](#2) Python library is used.

For the found structures, an ionic relaxation is performed for each of these structures with [M3GNet](https://github.com/materialsvirtuallab/m3gnet) [[3]](#3). M3GNet implement machine-learning interatomic potentials, thus relaxations can be performed in few seconds with a regular computer (instead of minutes or hours in a supercomputer). After this, it is strongly recommended to perform DFT relaxations with the most interesting resulting phases.

The relaxed structures are interesting enough to finish the process, but it is also possible to distort structures (with Gaussian noise) and relax them again in order to be able to reach unexplored low-energy phases. This is implemented with a simple Metropolis Algorithm, the ions of the less energetic phases after the initial relaxation are distorted and the new phases are relaxed. After these relaxations their energies are computed, and if are smaller than the original energies the new phase is accepted, if it is not we reject this phase. If we repeat this enough times (in what we call Generations Loop), thus, we have enough generations, we should be able to explore the unreached minimums in the energy hypersurface. 

It is also possible to look for crystal phases with a fixed pressure. The desired pressure can be introduced as an input and after the relaxation of the found structures, the program will determine the compressed structures and energies for the given pressure. For now, the program only takes into account compressions or expansions of the unit cell. Given a pressure the volume of the distorted structure can be found minimizing the enthalpy, $H$, that can be computed with:

$H(p,V)=E(V)+pV.$

After this we rank the pressurized structures as function of the enthalpy (instead of the internal energy as is the case of computations without pressure).

In order to perform pressure calculations you should use:

```
comp_pressure = True
```
Then, assign the desired values to the flags `pressure`, `num_volumes`, `minimum_volume` and `maximum_volume`.

Since it is possible to re-train M3GNet with DFT results, we can use DFT materials simulation results to re-train M3GNet in order to be able to look for novel phases of similar materials, but this is not yet implemented.

We can use the determined phases with PyMCSP to predict the phase group of a given experimental diffractogram. In order to do this, we should provide a CSV file with the intensities as function of $2\theta$, the diffraction script will compare how likely are the experimental curves $I_{exp}$ and the theoretical curves $I_{theo}$. In order to do this a loss factor is computed:

$Loss=a\left|n_{exp} - n_{theo} \right|^{2} + \int_{2\theta_{min}}^{2\theta_{max}}\left| I_{exp} -I_{theo} \right|^{2}d2\theta,$

here, the first term penalizes big difference of peak numbers, while the second term penalizes peaks in different $2\theta$ positions. After the execution of the program, the results are shown in a file. The program also slightly varies the lattice parameters, and takes the volume that minimize the loss factor. Doing this we can match the experimental lattice parameters.

## Requirements

The code has been tested in GNU/Linux operating system (Ubuntu 22.04.3 LTS) and with Python 3.10.12. The required packages to execute PyMCSP are: 
- numpy
- scipy
- pymatgen
- pyxtal
- m3gnet 0.2.4
- tensorflow 2.15.0.post1
- keras 2.15.0

The different modules can be downloaded manually, or instead execute the following command to install them automatically:
```bash
$ pip install -r requirements.txt
```

## Installation

To download the repository use:

```bash
$ git clone https://github.com/polbeni/PyMCSP
```

Make sure to have installed all the necessary packages.

## How to use it

Modify the values in the `src/inputs` file and execute the main script `src/PyMCSP.py` with:

```bash
$ cd src
$ python3 PyMCSP.py
```

The calculations will start and a directory named `structures_files` will storage the found structures as well as other output files. If you want to perform calculations with different inputs (for example, with different materials), you should save the `structures_files` directory in another place, since every time PyMCSP is executed it overwrites this directory.

To do diffraction studies, modify the values in the `src/inputs_diffraction` file, be sure you have a CSV file (with the proper format) with the diffraction results in the directory, and execute the script `src/diffraction.py` with:

```bash
$ cd src
$ python3 diffraction.py
```

It is also important that you previously generate structures with PyMCSP.

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


