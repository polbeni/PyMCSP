# PyMCSP

A Python and [M3GNet](https://github.com/materialsvirtuallab/m3gnet) [[1]](#1) implementation for Crystal Structure Prediction. It is intended to be a methodology that can find interesting and novel phases for materials in few computational time and in a personal computer.

The code is still in development, however it is functional and you can try it.

## How it works

From a given atoms and stoichiometry the program constructs random crystal structures with the phase groups that can be constructed with the given stoichiometry. To do this [PyXtal](https://github.com/qzhu2017/PyXtal) [[2]](#2) Python library is used.

After the generation of structures a relaxation is performed for each of them with the model [M3GNet](https://github.com/materialsvirtuallab/m3gnet). 

## Installation

To download the repository:

```bash
$ git clone https://github.com/polbeni/PyMCSP
```

## How to use it

Modify the values in the `src/inputs` file and execute the main script `src/PyMCSP.py` with:

```bash
$ cd src
$ python3 PyMCSP.py
```

## References
<a id="1">[1]</a> 
CHEN, Chi; ONG, Shyue Ping. A universal graph deep learning interatomic potential for the periodic table. <em>Nature Computational Science</em>, 2022, 2.11: 718-728.

<a id="2">[2]</a> 
FREDERICKS, Scott, et al. PyXtal: A Python library for crystal structure generation and symmetry analysis. <em>Computer Physics Communications</em>, 2021, 261: 107810.
