# Pol Benítez Colominas, date of creation: 2023/09/08, date of last modification: 2023/09/11

import os

from pyxtal import pyxtal

num_groups_0d = 58
num_groups_1d = 75
num_groups_2d = 80
num_groups_3d = 230

struc = pyxtal()

atoms = ['Ag', 'S', 'Br']
stoichiometry = [3, 1, 1]
num_max_atoms = 5

max_atoms_condition = False

stoi_counter = 1
stoichiometry_dict = {}
while max_atoms_condition == False:
    variable_name = f"stoichiometry_{stoi_counter}"

    stoichiometry_dict[variable_name] = [element * stoi_counter for element in stoichiometry]

    stoi_counter = stoi_counter + 1

    if sum(stoichiometry_dict[variable_name]) >= num_max_atoms:
        max_atoms_condition = True

os.mkdir('structure_files')
os.mkdir('structure_files/initial_structures')
os.mkdir('structure_files/initial_structures/generated_structures/')

num_phase = 1
for key, value in stoichiometry_dict.items():
    print(f"Stoichiometry: {value}")
    for num_space_group in range(num_groups_3d):
        try:
            struc.from_random(3, num_space_group+1, atoms, value)
            name_file = 'structure_files/initial_structures/generated_structures//POSCAR-' + "{:03d}".format(num_phase)
            struc.to_file(name_file, fmt='poscar')

            num_phase = num_phase + 1
        except Exception:
            print(f'No compatibility of the stoichiometry with the phase group {num_space_group + 1}')
            #pass
