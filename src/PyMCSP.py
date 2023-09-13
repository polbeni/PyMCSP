# Pol Benítez Colominas, Universitat Politècnica de Catalunya
# September 2023 - September 2023
# Version 0.2

# Main script to perfrom crystal structure prediction with PyMCSP method


#### Import libraries

import os
import shutil

import numpy as np

from functions_pymcsp import *
from terminal_outputs import *

from pymatgen.core import Structure
from pymatgen.io.vasp import Poscar

from pyxtal import pyxtal

import spglib

from m3gnet.models import Relaxer

import warnings

for category in (UserWarning, DeprecationWarning):
    warnings.filterwarnings("ignore", category=category, module="tensorflow")


#### Read inputs from inputs file

dimension = None
atoms = []
stoichiometry = []
num_max_atoms = None
num_same_phase = None
print_terminal_outputs = None
save_all_generations = None
structure_file = None
prec_group_det = None
surviving_phases = None
num_generations = None 
max_disp = None 

inputs = open('inputs', "r")

variables = [None]*12

for it in range(10):
    inputs.readline()

for it in range(len(variables)):
    line = inputs.readline()

    if line.split()[2] == 'True':
        variables[it] = True
    elif line.split()[2] == 'False':
        variables[it] = False
    elif (it == 1) or (it == 2):
        parts = line.strip().split('=')

        variables[it] = eval(parts[1].strip())
    else:
        variables[it] = line.split()[2]

inputs.close()

dimension = int(variables[0])
atoms = variables[1]
stoichiometry = variables[2]
num_max_atoms = int(variables[3])
num_same_phase = int(variables[4])
print_terminal_outputs = variables[5]
save_all_generations = variables[6]
structure_file = variables[7]
prec_group_det = float(variables[8])
surviving_phases = float(variables[9])
num_generations = int(variables[10])
max_disp = float(variables[11])


#### Generate the structures

max_atoms_condition = False

stoi_counter = 1
stoichiometry_dict = {}
while max_atoms_condition == False:
    variable_name = f"stoichiometry_{stoi_counter}"

    stoichiometry_dict[variable_name] = [element * stoi_counter for element in stoichiometry]

    stoi_counter = stoi_counter + 1

    if sum(stoichiometry_dict[variable_name]) >= num_max_atoms:
        max_atoms_condition = True

if os.path.exists('structure_files'):
    shutil.rmtree('structure_files')

if os.path.exists('final_structures'):
    shutil.rmtree('final_structures')

os.mkdir('structure_files')
os.mkdir('structure_files/initial_structures')
os.mkdir('structure_files/initial_structures/generated_structures/')

num_phase = 1
for num_phase_group in range(num_same_phase):
    total_num_phase = generate_phases(stoichiometry_dict, num_phase, dimension,
                                      atoms, print_terminal_outputs, structure_file)
    
    num_phase = total_num_phase


#### Relax the first structures

dir_path = 'structure_files/initial_structures/generated_structures/'

num_structures = 0

for path in os.listdir(dir_path):
    if os.path.isfile(os.path.join(dir_path, path)):
        num_structures = num_structures + 1

if print_terminal_outputs == True:
    print(f'Number of structures to relax: {num_structures}')

os.mkdir('structure_files/initial_structures/relaxed_structures/')

relaxer = Relaxer()

phase_energy_array = np.zeros((num_structures,2))
count = 0

for num_struc in range(num_structures):
    initial_path = 'structure_files/initial_structures/generated_structures/POSCAR-' + "{:04d}".format(num_struc + 1)
    relaxed_path = 'structure_files/initial_structures/relaxed_structures/POSCAR-' + "{:04d}".format(num_struc + 1)

    relax_energy, count = relax_structure(relaxer, initial_path, relaxed_path, 
                                          print_terminal_outputs, structure_file, count)
    
    phase_energy_array[count - 1,0] = int(count)
    phase_energy_array[count - 1,1] = relax_energy

sorted_indices = np.argsort(phase_energy_array[:,1])
phase_energy_sorted = phase_energy_array[sorted_indices]

file_energy = open('structure_files/initial_structures/relaxed_structures/energy_ranking.txt', 'w')
file_energy.write('#       POSCAR-num       energy per atom (eV)       Space Group (Hermann-Mauguin)\n')

count = 1
for num_struc in range(num_structures):
    struc_number = int(phase_energy_sorted[num_struc,0])
    struc_energy = phase_energy_sorted[num_struc,1]
    phase_file = 'structure_files/initial_structures/relaxed_structures/POSCAR-' + "{:04d}".format(struc_number)

    write_in_file(file_energy, phase_file, struc_number, struc_energy, 
                  prec_group_det, count)

    count = count + 1

file_energy.close()


#### Start the generation loop

files_in_dir = os.listdir('structure_files/initial_structures/relaxed_structures/')

number_ini_struc = 0
for filename in files_in_dir:
    if filename.startswith('POSCAR-'):
       number_ini_struc = number_ini_struc + 1

number_surv_struc = int(number_ini_struc*surviving_phases)

if num_generations != 0:

    for num_gen in range(num_generations):
        create_dir_generation(num_gen + 1)

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

            for num_str in range(number_surv_struc):
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

            relaxed_energy = relax_structure_gen(relaxer, path_dist_file, path_relax_file,
                                                print_terminal_outputs, structure_file)
            
            write_in_file_gen(file_energy, path_relax_file, struc_file, relaxed_energy, 
                            prec_group_det, count)
            
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

if (save_all_generations == True) and (num_generations != 0):
    os.mkdir('structure_files/final_structures')

    for struc_file in name_selected_structures:
        shutil.copy(f'structure_files/generation-{num_generations:03d}/initial_structures/' + struc_file, 
                    'structure_files/final_structures/' + struc_file)
    
    shutil.copy(f'structure_files/generation-{num_generations:03d}/final_structures/energy_ranking.txt', 
                    'structure_files/final_structures/energy_ranking.txt')
elif (save_all_generations == False) and (num_generations != 0):
    os.mkdir('final_structures')

    for struc_file in name_selected_structures:
        shutil.copy(f'structure_files/generation-{num_generations:03d}/initial_structures/' + struc_file, 
                    'final_structures/' + struc_file)
    
    shutil.copy(f'structure_files/generation-{num_generations:03d}/final_structures/energy_ranking.txt', 
                    'final_structures/energy_ranking.txt')
    
    shutil.rmtree('structure_files')