# PyMCSP

One of the most important and critical problems in materials science is to know the crystal structure or structures for a given material. The position that ions occupy in crystal determines some of the properties of this material, from if the material is an insulator, conductor or semiconductor to more exotic phenomena as superconductivity. Thus it is important to know the crystal structures of the different materials. Many experimental methodologies (as X-ray diffraction) help experimental scientists to determine these structures, but sometimes we could be interested in performing simulations of materials that have not been synthesized yet. Crystal structure prediction methods [[1]](#1), are methodologies that look for novel phases. These are methods that present very good performance, but that need to use supercomputer facilities and a lot of computational time, since they usually perform plenty of DFT simulations.

In this context, we present **PyMCSP** a **P**ython and [**M**3GNet](https://github.com/materialsvirtuallab/m3gnet) [[2]](#2) implementation for **C**rystal **S**tructure **P**rediction. This is intended to be a method to find novel phases for new (or not so new) materials in reasonable times (few minutes to few hours in a standard user laptop) and in a personal computer.

The code is still in development, however it is functional and you can try it.

## How it works

From a given atoms and stoichiometry, the program constructs random crystal structures with the phase groups that are compatible with the given stoichiometry. To do this [PyXtal](https://github.com/qzhu2017/PyXtal) [[3]](#3) Python library is used.

For the found structures, an ionic relaxation is performed for each of these structures with M3GNet. M3GNet implement machine-learning force fields, thus relaxations can be performed in few seconds with a regular computer (instead of minutes or hours in a supercomputer). The relaxed structures are interesting enough to finish the process, but it is also possible to distort structures and relax again in order to be able to reach less energetic phases (although this is still in development).

After this, it is strongly recommended to perform DFT relaxations with the most interesting phases.

Since it is possible to re-train M3GNet with DFT results, we can use DFT materials simulations results to re-train M3GNet in order to be able to look for novel phases of similar materials, but this is not yet implemented. 

Some examples of found structures can be find in `examples`. There are results for three different materials:
- Carbon (C)
- Ag<sub>3</sub>SBr
- BaTiO<sub>3</sub>

For each of these materials you can find the ten less energetic POSCAR files (by the energies provided by M3GNet), and the `inputs` file that we used to generate them.

## Requirments

The code has been tested in GNU/Linux operating system (Ubuntu 22.04.3 LTS) and with Python 3.10.12. The required packages to execute PyMCSP are: 
- numpy
- pymatgen
- pyxtal
- spglib
- m3gnet 0.2.4

## Installation

To download the repository:

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

The calculations will start and a directory named `structures_files` will storage the found structures as well as other output files. If you want to perfome calculations with different inputs (for example, with different materials), you should save the `structures_files` directory in another place, since everytime PyMCSP is executed it overwrites this directory.

With the actual version (0.2), it is strongly recommended to fix the value of the input variable `num_generations` in the file `inputs` as:

```
num_generations = 0
```

The reason is because for now the ion distort applies a gaussian noise and transforms the symmetry of the system in a triclinic cell. But in some situations we can be interested , thus the idea is to implement a method to distort the ions positions but mantaining the symmetry of the cell. We hope we can upgrade it soon.

Then, imposing `num_generations` equal to zero, we avoid larger computation times, since with the first simulation (generation of structures and M3GNet relaxation) we have enough to find some interesting phases.

## References

<a id="1">[1]</a> 
WOODLEY, Scott M.; CATLOW, Richard. Crystal structure prediction from first principles. <em>Nature materials</em>, 2008, 7.12: 937-946.

<a id="2">[2]</a> 
CHEN, Chi; ONG, Shyue Ping. A universal graph deep learning interatomic potential for the periodic table. <em>Nature Computational Science</em>, 2022, 2.11: 718-728.

<a id="3">[3]</a> 
FREDERICKS, Scott, et al. PyXtal: A Python library for crystal structure generation and symmetry analysis. <em>Computer Physics Communications</em>, 2021, 261: 107810.
