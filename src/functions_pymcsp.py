# Pol Benítez Colominas, Universitat Politècnica de Catalunya
# September 2023 - September 2023
# Version 0.2

# Functions file

import os

import numpy as np

from pymatgen.core import Structure
from pymatgen.io.vasp import Poscar

from pyxtal import pyxtal

import spglib

def generate_phases(stoi_dict, phase_number, dim, atoms_arr, text_output, type_output):
    """
    Generate the initial phases

    Inputs:
        stoi_dict: the dictionary with the different stoichiometries that verifies
                   that number of atoms in unit cell <= num_max_atoms
        phase_number: name of the actual phase
        dim: dimension of the system (2 or 3)
        atoms_arr: array with atoms
        text_output: to decide if we want output text in terminal or not
        type_output: the crystal type of output file
    """

    num_groups = {
        2: 80,
        3: 230
    }

    struc = pyxtal()

    for key, value in stoi_dict.items():
        if text_output == True:
            print(f"Stoichiometry: {value}")

        for num_space_group in range(num_groups[dim]):
            try:
                struc.from_random(dim, num_space_group+1, atoms_arr, value)
                name_file = 'structure_files/initial_structures/generated_structures/POSCAR-' + "{:04d}".format(phase_number)
                struc.to_file(name_file, fmt=type_output)

                phase_number = phase_number + 1
            except Exception:
                if text_output == True:
                    print(f'No compatibility of the stoichiometry with the phase group {num_space_group + 1}')
                else:
                    pass

    return phase_number


def relax_structure(relax_object, init_struc_path, relax_struc_path, text_output, type_output, counter):
    """
    Relax an structure and get the energy

    Inputs:
        relax_object: the relaxer object by m3gnet module
        init_struc_path: the path of the structure to relax
        relax_struc_path: the path of the relaxed structure
        text_output: to decide if we want output text in terminal or not
        type_output: the crystal type of output file
        counter: a number counter for the phase
    """

    crystal_structure = Structure.from_file(init_struc_path)

    if text_output == True:
        relaxed_structure = relax_object.relax(crystal_structure, verbose=True)
    else:
        relaxed_structure = relax_object.relax(crystal_structure, verbose=False)

    final_structure = relaxed_structure['final_structure']

    final_structure.to(filename=relax_struc_path, fmt=type_output)

    final_energy = float(relaxed_structure['trajectory'].energies[-1]/len(crystal_structure))
    counter = counter + 1

    return final_energy, counter 


def write_in_file(file_object, structure_path, struc_num, energy, prec_spglib, counter):
    """
    Write the number, the file structure number, the energy per atom (in eV),
    and the phase group for a given structure in the energy file

    Inputs:
        file_object: the object of the file of energies
        structure_path: path to the structure of interest
        struc_num: the number of the structure
        energy: energy per atom (in eV) of the structure
        prec_spglib: precision for the phase determination with spglib
        counter: a number counter for the phase
    """

    structure = Poscar.from_file(structure_path).structure

    lattice = structure.lattice.matrix
    positions = [site.coords for site in structure.sites]
    atomic_numbers = [site.specie.number for site in structure.sites]

    spglib_cell = (lattice, positions, atomic_numbers)

    space_group_info = spglib.get_spacegroup(spglib_cell, symprec=prec_spglib, 
                                             angle_tolerance=-1.0, symbol_type=0)
    
    file_object.write(f'{int(counter)}       POSCAR-{struc_num:04d}       {energy}       {space_group_info}\n')

    return


def create_dir_generation(number_generation):
    """
    Creates the necessary directories to store all the files for each generation

    Inputs:
        number_generation: number of the actual generation
    """

    os.mkdir(f'structure_files/generation-{number_generation:03d}')
    os.mkdir(f'structure_files/generation-{number_generation:03d}/initial_structures')
    os.mkdir(f'structure_files/generation-{number_generation:03d}/distorsed_structures')
    os.mkdir(f'structure_files/generation-{number_generation:03d}/relaxed_structures')
    os.mkdir(f'structure_files/generation-{number_generation:03d}/final_structures')

    return


def structure_distortion(file, max_displacement, final_path):
    """
    Distort a crystal structure

    Inputs:
        file: name or path of the file
        max_displacement: maximum displacement in a given direction (in Angstroms)
        final_path: name and path of the final file
    """

    structure = Structure.from_file(file)

    displaced_structure = structure.copy()
    for site in displaced_structure:
        displacement_vector = np.random.uniform(-max_displacement, max_displacement, 3)
        site.coords = site.coords + displacement_vector

    displaced_structure.to(filename=final_path, fmt='poscar')

    return


def relax_structure_gen(relax_object, dist_struc_path, relax_struc_path, text_output, type_output):
    """
    Relax an structure and get the energy for the generations part of the main code

    Inputs:
        relax_object: the relaxer object by m3gnet module
        dist_struc_path: the path of the distorted structure
        relax_struc_path: the path of the relaxed structure
        text_output: to decide if we want output text in terminal or not
        type_output: the crystal type of output file
    """

    crystal_structure = Structure.from_file(dist_struc_path)

    if text_output == True:
        relaxed_structure = relax_object.relax(crystal_structure, verbose=True)
    else:
        relaxed_structure = relax_object.relax(crystal_structure, verbose=False)

    final_structure = relaxed_structure['final_structure']

    final_structure.to(filename=relax_struc_path, fmt=type_output)

    final_energy = float(relaxed_structure['trajectory'].energies[-1]/len(crystal_structure))


    return final_energy


def write_in_file_gen(file_object, structure_path, struc_num, energy, prec_spglib, counter):
    """
    Write the number, the file structure number, the energy per atom (in eV),
    and the phase group for a given structure in the energy file

    Inputs:
        file_object: the object of the file of energies
        structure_path: path to the structure of interest
        struc_num: the number of the structure
        energy: energy per atom (in eV) of the structure
        prec_spglib: precision for the phase determination with spglib
        counter: a number counter for the phase
    """

    structure = Poscar.from_file(structure_path).structure

    lattice = structure.lattice.matrix
    positions = [site.coords for site in structure.sites]
    atomic_numbers = [site.specie.number for site in structure.sites]

    spglib_cell = (lattice, positions, atomic_numbers)

    space_group_info = spglib.get_spacegroup(spglib_cell, symprec=prec_spglib, 
                                             angle_tolerance=-1.0, symbol_type=0)
    
    file_object.write(f'{int(counter)}       {struc_num}       {energy}       {space_group_info}\n')

    return