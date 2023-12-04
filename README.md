# PyMCSP

One of the most important and critical problems in materials science is to know the crystal structure or structures for a given material. The position that ions occupy in the unit cell of a crystal and in the periodic table of elements fully determines the physical, chemical and functional properties of the material, from if the material is an insulator, conductor or semiconductor to more exotic phenomena as superconductivity. Thus, it is important to know the crystal structures of the different materials. Scattering experiments (as X-ray diffraction) help experimental scientists to determine these structures, but sometimes we could be interested in performing simulations of materials that have not been synthesized yet and of which we know nothing of their structural properties. Crystal structure prediction methods [[1]](#1), are methodologies that look for local (and hopefully the global) minimums of energy, thus for the phases. These methods usually present a good performance and in some cases are able to find novel phases, but that need to use supercomputer facilities and a lot of computer time, since they usually perform plenty of DFT simulations.

In this context, we present **PyMCSP** a **P**ython and **M**achine Learning methods implementation for **C**rystal **S**tructure **P**rediction. This is intended to be a method to find novel phases for well known or new materials at reasonable times, few minutes to few hours in a standard user's laptop.

The code is still in development and many new implementations are intended to be added, however it is functional and you can try it.

## How it works

From a given chemical compounds and stoichiometry, the program constructs random crystal structures with the phase groups that are compatible with the given stoichiometry. To do this [PyXtal](https://github.com/qzhu2017/PyXtal) [[2]](#2) Python library is used.

For the found structures, an ionic relaxation is performed for each of these structures with [M3GNet](https://github.com/materialsvirtuallab/m3gnet) [[3]](#3). M3GNet implement machine-learning interatomic potentials, thus relaxations can be performed in few seconds with a regular computer (instead of minutes or hours in a supercomputer). The relaxed structures are interesting enough to finish the process, but it is also possible to distort structures and relax again in order to be able to reach less energetic phases (although, this is still in development and some problems can arise if it is used).

After this, it is strongly recommended to perform DFT relaxations with the most interesting resulting phases.

Since it is possible to re-train M3GNet with DFT results, we can use DFT materials simulation results to re-train M3GNet in order to be able to look for novel phases of similar materials, but this is not yet implemented (expected to be implemented soon). 

Some examples of resulting structures can be found in `examples` directory. There are results for three different materials:
- Carbon (C)
- Ag<sub>3</sub>SBr
- BN

For each of these materials you can find the ten less energetic POSCAR files (by the energies provided by M3GNet), and the `inputs` file that we used to generate them.

## Requirments

The code has been tested in GNU/Linux operating system (Ubuntu 22.04.3 LTS) and with Python 3.10.12. The required packages to execute PyMCSP are: 
- numpy
- pymatgen
- pyxtal
- m3gnet 0.2.4

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

With the actual version (0.2), it is strongly recommended to fix the value of the input variable `num_generations` in the file `inputs` as:

```
num_generations = 0
```

The reason of this is because for now the ion distortion applies a Gaussian noise and transforms the symmetry of the system in a triclinic cell. 

Then, imposing `num_generations` equal to zero, we avoid larger computation times, since with the first simulation (generation of structures and M3GNet relaxation) we have enough to find some interesting phases.

## References

<a id="1">[1]</a> 
WOODLEY, Scott M.; CATLOW, Richard. Crystal structure prediction from first principles. <em>Nature materials</em>, 2008, 7.12: 937-946.

<a id="2">[2]</a> 
FREDERICKS, Scott, et al. PyXtal: A Python library for crystal structure generation and symmetry analysis. <em>Computer Physics Communications</em>, 2021, 261: 107810.

<a id="3">[3]</a> 
CHEN, Chi; ONG, Shyue Ping. A universal graph deep learning interatomic potential for the periodic table. <em>Nature Computational Science</em>, 2022, 2.11: 718-728.


