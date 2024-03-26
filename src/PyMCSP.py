# Pol Benítez Colominas, Universitat Politècnica de Catalunya
# September 2023 - March 2024
# Version 0.3

# Main script to perfrom crystal structure prediction with PyMCSP method


#### Import libraries

import os
import sys
import shutil
import time

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

star_time = time.time()

#### Read inputs from inputs file

dimension = None
atoms = []
stoichiometry = []
num_max_atoms = None
num_same_phase = None
max_ionic_steps = None
comp_pressure = None
pressure = None
num_volumes = None
minimum_volume = None 
maximum_volume = None
print_terminal_outputs = None
save_log_file = None
save_all_generations = None
structure_file = None
prec_group_det = None
num_generations = None
surviving_phases = None
max_disp = None 

inputs = open('inputs', "r")

variables = [None]*19

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
max_ionic_steps = int(variables[5])
comp_pressure = variables[6]
pressure = float(variables[7])
num_volumes = int(variables[8])
minimum_volume = float(variables[9]) 
maximum_volume = float(variables[10])
print_terminal_outputs = variables[11]
save_log_file = variables[12]
save_all_generations = variables[13]
structure_file = variables[14]
prec_group_det = float(variables[15])
num_generations = int(variables[16])
surviving_phases = float(variables[17])
max_disp = float(variables[18])


### Save the log file
    
log_file = 'output.log'

if save_log_file == True:
    log_file_handle = open(log_file, "w")

    class Tee:
        def __init__(self, *files):
            self.files = files

        def write(self, text):
            for file in self.files:
                file.write(text)
                file.flush()

        def flush(self):
            for file in self.files:
                file.flush()

    sys.stdout = Tee(sys.stdout, log_file_handle)

if print_terminal_outputs == True:
    initial_message()

    simulation_information(atoms, stoichiometry, num_max_atoms, num_generations)


#### Generate the structures

if print_terminal_outputs == True:
    struc_gen_ini()

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

if print_terminal_outputs == True:
    struc_gen_end()


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

if print_terminal_outputs == True:
    relax_ini()

for num_struc in range(num_structures):
    if structure_file == 'poscar':
        initial_path = 'structure_files/initial_structures/generated_structures/POSCAR-' + "{:06d}".format(num_struc + 1)
        relaxed_path = 'structure_files/initial_structures/relaxed_structures/POSCAR-' + "{:06d}".format(num_struc + 1)
    elif structure_file == 'cif':
        initial_path = 'structure_files/initial_structures/generated_structures/structure-' + "{:06d}".format(num_struc + 1) + '.cif'
        relaxed_path = 'structure_files/initial_structures/relaxed_structures/structure-' + "{:06d}".format(num_struc + 1) + '.cif'

    relax_energy, count = relax_structure(relaxer, initial_path, relaxed_path, max_ionic_steps, 
                                          print_terminal_outputs, structure_file, count)
    
    phase_energy_array[count - 1, 0] = int(count)
    phase_energy_array[count - 1, 1] = relax_energy

if print_terminal_outputs == True:
    relax_end()

sorted_indices = np.argsort(phase_energy_array[:,1])
phase_energy_sorted = phase_energy_array[sorted_indices]

file_energy = open('structure_files/initial_structures/relaxed_structures/energy_ranking.txt', 'w')
if structure_file == 'poscar':
    file_energy.write('#       POSCAR-num       energy per atom (eV)       Space Group (Hermann-Mauguin)\n')
elif structure_file == 'cif':
    file_energy.write('#       structure-num.cif       energy per atom (eV)       Space Group (Hermann-Mauguin)\n')

count = 1
for num_struc in range(num_structures):
    struc_number = int(phase_energy_sorted[num_struc,0])
    struc_energy = phase_energy_sorted[num_struc,1]
    if structure_file == 'poscar':
        phase_file = 'structure_files/initial_structures/relaxed_structures/POSCAR-' + "{:06d}".format(struc_number)
    elif structure_file == 'cif':
        phase_file = 'structure_files/initial_structures/relaxed_structures/structure-' + "{:06d}".format(struc_number) + '.cif'

    write_in_file(file_energy, phase_file, struc_number, struc_energy, 
                  prec_group_det, structure_file, count)

    count = count + 1

file_energy.close()


### Pressure computations

if comp_pressure == True:
    dir_path = 'structure_files/initial_structures/relaxed_structures/'

    num_structures = 0

    for path in os.listdir(dir_path):
        if os.path.isfile(os.path.join(dir_path, path)):
            num_structures = num_structures + 1

    if print_terminal_outputs == True:
        print(f'Number of structures to study with pressure conditions: {num_structures - 1}')

    os.mkdir('structure_files/initial_structures/pressure_structures/')

    relaxer = Relaxer()

    pressure_energy_array = np.zeros((num_structures - 1,3))
    count = 0

    if print_terminal_outputs == True:
        pressure_ini()

    for num_struc in range(num_structures - 1):
        if structure_file == 'poscar':
            relaxed_path = 'structure_files/initial_structures/relaxed_structures/POSCAR-' + "{:06d}".format(num_struc + 1)
            pressure_path = 'structure_files/initial_structures/pressure_structures/POSCAR-' + "{:06d}".format(num_struc + 1)
        elif structure_file == 'cif':
            relaxed_path = 'structure_files/initial_structures/relaxed_structures/structure-' + "{:06d}".format(num_struc + 1) + '.cif'
            pressure_path = 'structure_files/initial_structures/pressure_structures/structure-' + "{:06d}".format(num_struc + 1) + '.cif'

        pressure_energy, count, alpha = pressure_structure(relaxer, relaxed_path, pressure_path, pressure, num_volumes, minimum_volume,
                                                           maximum_volume, print_terminal_outputs, structure_file, count)
        
        pressure_energy_array[count - 1, 0] = int(count)
        pressure_energy_array[count - 1, 1] = pressure_energy
        pressure_energy_array[count - 1, 2] = alpha

    if print_terminal_outputs == True:
        pressure_end()

    sorted_indices_pressure = np.argsort(pressure_energy_array[:,1])
    phase_energy_sorted_pressure = pressure_energy_array[sorted_indices_pressure]

    file_energy = open('structure_files/initial_structures/pressure_structures/energy_ranking_pressure.txt', 'w')
    if structure_file == 'poscar':
        file_energy.write('#       POSCAR-num       energy per atom (eV)       alpha       Space Group (Hermann-Mauguin)\n')
    elif structure_file == 'cif':
        file_energy.write('#       structure-num.cif       energy per atom (eV)       alpha       Space Group (Hermann-Mauguin)\n')

    count = 1
    for num_struc in range(num_structures - 1):
        struc_number = int(phase_energy_sorted_pressure[num_struc,0])
        struc_energy = phase_energy_sorted_pressure[num_struc,1]
        struc_alpha = phase_energy_sorted_pressure[num_struc,2]
        if structure_file == 'poscar':
            phase_file = 'structure_files/initial_structures/pressure_structures/POSCAR-' + "{:06d}".format(struc_number)
        elif structure_file == 'cif':
            phase_file = 'structure_files/initial_structures/pressure_structures/structure-' + "{:06d}".format(struc_number) + '.cif'

        if struc_energy == 0:
            struc_energy = 'Error'
            struc_alpha = 'Not minimized'

        write_in_file_pressure(file_energy, phase_file, struc_number, struc_energy, struc_alpha,
                    prec_group_det, structure_file, count)

        count = count + 1

    file_energy.close()


#### Start the generation loop

files_in_dir = os.listdir('structure_files/initial_structures/relaxed_structures/')

number_ini_struc = 0
if structure_file == 'poscar':
    for filename in files_in_dir:
        if filename.startswith('POSCAR-'):
            number_ini_struc = number_ini_struc + 1
elif structure_file == 'cif':
    for filename in files_in_dir:
        if filename.startswith('structure-'):
            number_ini_struc = number_ini_struc + 1

number_surv_struc = int(number_ini_struc*surviving_phases)

if num_generations != 0:

    if print_terminal_outputs == True:
        gen_ini()

    for num_gen in range(num_generations):

        if print_terminal_outputs == True:
            gen_ini_actual(num_gen + 1)

        create_dir_generation(num_gen + 1)

        name_dir_initial_struc = {
            'initial': 'structure_files/initial_structures/relaxed_structures/',
            'not_initial': f'structure_files/generation-{num_gen:03d}/final_structures'
        }

        path_destination = f'structure_files/generation-{num_gen + 1:03d}/initial_structures/'
        if structure_file == 'poscar':
            file_prefix = 'POSCAR-'
        elif structure_file == 'cif':
            file_prefix = 'structure-'

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
            path_dist_file = f'structure_files/generation-{num_gen + 1:03d}/distorted_structures/' + struc_file

            structure_distortion(path_struc_file, max_disp, structure_file, path_dist_file)

        file_energy = open(f'structure_files/generation-{num_gen + 1:03d}/relaxed_structures/energy_ranking.txt', 'w')
        if structure_file == 'poscar':
            file_energy.write('#       POSCAR-num       energy per atom (eV)       Space Group (Hermann-Mauguin)\n')
        elif structure_file == 'cif':
            file_energy.write('#       structure-num.cif       energy per atom (eV)       Space Group (Hermann-Mauguin)\n')
        count = 1

        if print_terminal_outputs == True:
            relax_ini()

        for struc_file in name_selected_structures:
            path_dist_file = f'structure_files/generation-{num_gen + 1:03d}/distorted_structures/' + struc_file
            path_relax_file = f'structure_files/generation-{num_gen + 1:03d}/relaxed_structures/' + struc_file

            relaxed_energy = relax_structure_gen(relaxer, path_dist_file, path_relax_file, max_ionic_steps,
                                                print_terminal_outputs, structure_file)
            
            write_in_file_gen(file_energy, path_relax_file, struc_file, relaxed_energy, 
                            prec_group_det, structure_file, count)
            
            count = count + 1

        if print_terminal_outputs == True:
            relax_end()

        file_energy.close()

        path_previous_energy_file = {
            'initial': 'structure_files/initial_structures/relaxed_structures/energy_ranking.txt',
            'not_initial': f'structure_files/generation-{num_gen:03d}/final_structures/energy_ranking.txt'
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
        if structure_file == 'poscar':
            file_energy.write('#       POSCAR-num       energy per atom (eV)       Space Group (Hermann-Mauguin)\n')
        elif structure_file == 'cif':
            file_energy.write('#       structure-num.cif       energy per atom (eV)       Space Group (Hermann-Mauguin)\n')

        for struc_file in name_selected_structures:
            line_previous_energy = previous_energy.readline()
            line_actual_energy = actual_energy.readline()

            struc_previous_energy = float(line_previous_energy.split()[2])
            struc_actual_energy = float(line_actual_energy.split()[2])

            if (struc_previous_energy >= struc_actual_energy) and (abs(struc_previous_energy - struc_actual_energy) >= 2e-4):
                shutil.copy(f'structure_files/generation-{num_gen + 1:03d}/relaxed_structures/' + struc_file, 
                            f'structure_files/generation-{num_gen + 1:03d}/final_structures')
                
                file_energy.write(line_actual_energy)
            else:
                shutil.copy(f'structure_files/generation-{num_gen + 1:03d}/initial_structures/' + struc_file, 
                            f'structure_files/generation-{num_gen + 1:03d}/final_structures')
                
                file_energy.write(line_previous_energy)

        file_energy.close()        
        previous_energy.close()
        actual_energy.close()

        if print_terminal_outputs == True:
            gen_end_actual(num_gen + 1)

        num_gen = num_gen + 1

    if print_terminal_outputs == True:
        gen_end()

if (save_all_generations == True) and (num_generations != 0):
    os.mkdir('structure_files/final_structures')

    for struc_file in name_selected_structures:
        shutil.copy(f'structure_files/generation-{num_generations:03d}/initial_structures/' + struc_file, 
                    'structure_files/final_structures/' + struc_file)
        
    phases_and_energies = []
    final_energy_file = open(f'structure_files/generation-{num_generations:03d}/final_structures/energy_ranking.txt', 'r')
    final_energy_file.readline()
    for it in range(len(name_selected_structures)):
        line_file = final_energy_file.readline()
        phases_and_energies_element = [line_file.split()[1], float(line_file.split()[2]), line_file.split()[3] + ' ' + line_file.split()[4]]
        phases_and_energies.append(phases_and_energies_element)
    final_energy_file.close()
    
    sorted_phases_and_energies = sorted(phases_and_energies, key=lambda x: x[1])

    final_energy_file = open('structure_files/final_structures/energy_ranking.txt', 'w')
    if structure_file == 'poscar':
        final_energy_file.write('#       POSCAR-num       energy per atom (eV)       Space Group (Hermann-Mauguin)\n')
    elif structure_file == 'cif':
        final_energy_file.write('#       structure-num.cif       energy per atom (eV)       Space Group (Hermann-Mauguin)\n')

    num_phase = 1
    for phase in sorted_phases_and_energies:
        final_energy_file.write(f'{num_phase}       {phase[0]}       {phase[1]}       {phase[2]}\n')
        num_phase = num_phase + 1
    
    final_energy_file.close()
elif (save_all_generations == False) and (num_generations != 0):
    os.mkdir('final_structures')

    for struc_file in name_selected_structures:
        shutil.copy(f'structure_files/generation-{num_generations:03d}/initial_structures/' + struc_file, 
                    'final_structures/' + struc_file)
        
    phases_and_energies = []
    final_energy_file = open(f'structure_files/generation-{num_generations:03d}/final_structures/energy_ranking.txt', 'r')
    final_energy_file.readline()
    for it in range(len(name_selected_structures)):
        line_file = final_energy_file.readline()
        phases_and_energies_element = [line_file.split()[1], float(line_file.split()[2]), line_file.split()[3] + ' ' + line_file.split()[4]]
        phases_and_energies.append(phases_and_energies_element)
    final_energy_file.close()
    
    sorted_phases_and_energies = sorted(phases_and_energies, key=lambda x: x[1])

    final_energy_file = open('final_structures/energy_ranking.txt', 'w')
    if structure_file == 'poscar':
        final_energy_file.write('#       POSCAR-num       energy per atom (eV)       Space Group (Hermann-Mauguin)\n')
    elif structure_file == 'cif':
        final_energy_file.write('#       structure-num.cif       energy per atom (eV)       Space Group (Hermann-Mauguin)\n')

    num_phase = 1
    for phase in sorted_phases_and_energies:
        final_energy_file.write(f'{num_phase}       {phase[0]}       {phase[1]}       {phase[2]}\n')
        num_phase = num_phase + 1
    
    final_energy_file.close()
    
    shutil.rmtree('structure_files')

end_time = time.time()
elapsed_time = end_time - star_time

if print_terminal_outputs == True:
    print(f'Elapsed time: {elapsed_time:.3f} seconds')

if print_terminal_outputs == True:
            final_message()


### Close the log file
    
if save_log_file == True:
    log_file_handle.close()
    sys.stdout = sys.__stdout__