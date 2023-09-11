# Pol Benítez Colominas, date of creation: 2023/09/11, date of last modification: 2023/09/11

import os
import shutil

import numpy as np

from pymatgen.core import Structure

import warnings

from m3gnet.models import Relaxer
from pymatgen.core import Lattice, Structure

from pymatgen.io.vasp import Poscar
import spglib

def structure_distortion(file, max_displacement, final_path):
    """
    Distort a crystal structure

    file: name or path of the file
    max_displacement: maximum displacement in a given direction (in Angstroms)
    final_path: name and path of the final file
    """

    structure = Structure.from_file(file)

    displaced_structure = structure.copy()  # Create a copy to avoid modifying the original
    for site in displaced_structure:
        displacement_vector = np.random.uniform(-max_displacement, max_displacement, 3)
        site.coords = site.coords + displacement_vector

    displaced_structure.to(filename=final_path, fmt='poscar')

    return

surviving_phases = 0.15
num_generations = 4
max_disp = 0.1

files_in_dir = os.listdir('structure_files/initial_structures/relaxed_structures/')

number_ini_struc = 0
for filename in files_in_dir:
    if filename.startswith('POSCAR-'):
       number_ini_struc = number_ini_struc + 1


number_struc = int(number_ini_struc*surviving_phases)

relaxer = Relaxer()

for num_gen in range(num_generations):
    if os.path.exists(f'structure_files/generation-{num_gen + 1:03d}'):
        shutil.rmtree(f'structure_files/generation-{num_gen + 1:03d}')


    os.mkdir(f'structure_files/generation-{num_gen + 1:03d}')
    os.mkdir(f'structure_files/generation-{num_gen + 1:03d}/initial_structures')
    os.mkdir(f'structure_files/generation-{num_gen + 1:03d}/distorsed_structures')
    os.mkdir(f'structure_files/generation-{num_gen + 1:03d}/relaxed_structures')
    os.mkdir(f'structure_files/generation-{num_gen + 1:03d}/final_structures')

    name_dir_initial_struc = {
        'initial': 'structure_files/initial_structures/relaxed_structures/',
        'not_initial': f'structure_files/generation-{num_gen:03d}/final_structures'
    }

    path_destination = f'structure_files/generation-{num_gen + 1:03d}/initial_structures/'

    file_prefix = 'POSCAR-'

    if  num_gen == 0:
        path_initial_struc = name_dir_initial_struc['initial']

        name_selected_structures = []

        energy_file = open('structure_files/initial_structures/relaxed_structures/energy_ranking.txt', "r")

        energy_file.readline()

        for num_str in range(number_struc):
            actual_line = energy_file.readline()

            name_selected_structures.append(actual_line.split()[1])
        

        energy_file.close()
        
        for filename in name_selected_structures:
            if filename.startswith(file_prefix):
                source_file = os.path.join(path_initial_struc, filename)
                destination_file = os.path.join(path_destination, filename)
                shutil.copy(source_file, destination_file)

    else:
        path_initial_struc = name_dir_initial_struc['not_initial']

        files_in_dir = os.listdir(path_initial_struc)

        for filename in files_in_dir:
            if filename.startswith(file_prefix):
                source_file = os.path.join(path_initial_struc, filename)
                destination_file = os.path.join(path_destination, filename)
                shutil.copy(source_file, destination_file)



    for struc_file in name_selected_structures:
        path_struc_file = f'structure_files/generation-{num_gen + 1:03d}/initial_structures/' + struc_file
        path_dist_file = f'structure_files/generation-{num_gen + 1:03d}/distorsed_structures/' + struc_file

        structure_distortion(path_struc_file, max_disp, path_dist_file)

    file_energy = open(f'structure_files/generation-{num_gen + 1:03d}/relaxed_structures/energy_ranking.txt', 'w')
    file_energy.write('#       POSCAR-num       energy per atom (eV)       Space Group (Hermann-Mauguin)\n')
    count = 1
    for struc_file in name_selected_structures:
        path_dist_file = f'structure_files/generation-{num_gen + 1:03d}/distorsed_structures/' + struc_file
        path_relax_file = f'structure_files/generation-{num_gen + 1:03d}/relaxed_structures/' + struc_file

        crystal_structure = Structure.from_file(path_dist_file)

        relaxed_structure = relaxer.relax(crystal_structure, verbose=True)

        final_structure = relaxed_structure['final_structure']

        final_structure.to(filename=path_relax_file, fmt='poscar')

        relaxed_energy = float(relaxed_structure['trajectory'].energies[-1]/len(crystal_structure))

        lattice = final_structure.lattice.matrix
        positions = [site.coords for site in final_structure.sites]
        atomic_numbers = [site.specie.number for site in final_structure.sites]

        spglib_cell = (lattice, positions, atomic_numbers)

        space_group_info = spglib.get_spacegroup(spglib_cell, symprec=1e-5, angle_tolerance=-1.0, symbol_type=0)


        file_energy.write(f'{int(count)}       {struc_file}       {relaxed_energy}       {space_group_info}\n')

        count = count + 1

    file_energy.close()

    path_previous_energy_file = {
        'initial': 'structure_files/initial_structures/relaxed_structures/energy_ranking.txt',
        'not_initial': f'structure_files/generation-{num_gen:03d}/relaxed_structures/energy_ranking.txt'
    }

    if num_gen == 0:
        previous_energy_path = path_previous_energy_file['initial']
    else:
        previous_energy_path = path_previous_energy_file['not_initial']

    actual_energy_path = f'structure_files/generation-{num_gen + 1:03d}/relaxed_structures/energy_ranking.txt'

    previous_energy = open(previous_energy_path, "r")
    previous_energy.readline()

    actual_energy = open(actual_energy_path, "r")
    actual_energy.readline()

    file_energy = open(f'structure_files/generation-{num_gen + 1:03d}/final_structures/energy_ranking.txt', 'w')
    file_energy.write('#       POSCAR-num       energy per atom (eV)       Space Group (Hermann-Mauguin)\n')

    for struc_file in name_selected_structures:
        line_previous_energy = previous_energy.readline()
        line_actual_energy = actual_energy.readline()

        struc_previous_energy = float(line_previous_energy.split()[2])
        struc_actual_energy = float(line_actual_energy.split()[2])

        if struc_previous_energy <= struc_actual_energy:
            shutil.copy(f'structure_files/generation-{num_gen + 1:03d}/initial_structures/' + struc_file, 
                        f'structure_files/generation-{num_gen + 1:03d}/final_structures')
            
            file_energy.write(line_previous_energy)
        else:
            shutil.copy(f'structure_files/generation-{num_gen + 1:03d}/relaxed_structures/' + struc_file, 
                        f'structure_files/generation-{num_gen + 1:03d}/final_structures')
            
            file_energy.write(line_actual_energy)
            
    file_energy.close()        
    previous_energy.close()
    actual_energy.close()
    

    num_gen = num_gen + 1

if os.path.exists('structure_files/final_structures'):
    shutil.rmtree('structure_files/final_structures')

os.mkdir('structure_files/final_structures')

for struc_file in name_selected_structures:
    shutil.copy(f'structure_files/generation-{num_generations:03d}/initial_structures/' + struc_file, 
                'structure_files/final_structures/' + struc_file)
shutil.copy(f'structure_files/generation-{num_generations:03d}/final_structures/energy_ranking.txt', 
                'structure_files/final_structures/energy_ranking.txt')