# Pol Benítez Colominas, date of creation: 2023/09/08, date of last modification: 2023/09/11

import os

import numpy as np

from pymatgen.core import Structure

import warnings

from m3gnet.models import Relaxer
from pymatgen.core import Lattice, Structure

from pymatgen.io.vasp import Poscar
import spglib

for category in (UserWarning, DeprecationWarning):
    warnings.filterwarnings("ignore", category=category, module="tensorflow")

dir_path = 'structure_files/initial_structures/generated_structures/'

num_POSCARs = 0

for path in os.listdir(dir_path):
    if os.path.isfile(os.path.join(dir_path, path)):
        num_POSCARs = num_POSCARs + 1

print(f'Number of POSCAR files to relax: {num_POSCARs}')

os.mkdir('structure_files/initial_structures/relaxed_structures/')

relaxer = Relaxer()

phase_energy_array = np.zeros((num_POSCARs,2))
count = 0

for num_POSCAR in range(num_POSCARs):
    POSCAR_file = 'structure_files/initial_structures/generated_structures/POSCAR-' + "{:03d}".format(num_POSCAR + 1)
    CONTCAR_file = 'structure_files/initial_structures/relaxed_structures/POSCAR-' + "{:03d}".format(num_POSCAR + 1)

    crystal_structure = Structure.from_file(POSCAR_file)

    relaxed_POSCAR = relaxer.relax(crystal_structure, verbose=True)

    final_POSCAR = relaxed_POSCAR['final_structure']

    final_POSCAR.to(filename=CONTCAR_file, fmt='poscar')

    phase_energy_array[count,0] = int(count + 1)
    phase_energy_array[count,1] = float(relaxed_POSCAR['trajectory'].energies[-1]/len(crystal_structure))
    
    count = count + 1

sorted_indices = np.argsort(phase_energy_array[:,1])
phase_energy_sorted = phase_energy_array[sorted_indices]

file_energy = open('structure_files/initial_structures/relaxed_structures/energy_ranking.txt', 'w')
file_energy.write('#       POSCAR-num       energy per atom (eV)       Space Group (Hermann-Mauguin)\n')
count = 1
for x in range(num_POSCARs): 
    phase_file = 'structure_files/initial_structures/relaxed_structures/POSCAR-' + "{:03d}".format(int(phase_energy_sorted[x,0]))
    structure = Poscar.from_file(phase_file).structure

    lattice = structure.lattice.matrix
    positions = [site.coords for site in structure.sites]
    atomic_numbers = [site.specie.number for site in structure.sites]

    spglib_cell = (lattice, positions, atomic_numbers)

    space_group_info = spglib.get_spacegroup(spglib_cell, symprec=1e-5, angle_tolerance=-1.0, symbol_type=0)

    file_energy.write(f'{int(count)}       POSCAR-{int(phase_energy_sorted[x,0]):03d}       {phase_energy_sorted[x,1]}       {space_group_info}\n')
        
    count = count + 1
file_energy.close()