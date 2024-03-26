# Pol Benítez Colominas, Universitat Politècnica de Catalunya
# September 2023 - March 2024
# Version 0.3

# Functions file

import os
from copy import deepcopy

import numpy as np
from scipy.interpolate import interp1d

from pymatgen.core import Structure
from pymatgen.io.vasp import Poscar
from pymatgen.io.cif import CifParser
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from pyxtal import pyxtal


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
        type_output: the crystal type of output file, poscar or cif
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
                if type_output == 'poscar':
                    name_file = 'structure_files/initial_structures/generated_structures/POSCAR-' + "{:06d}".format(phase_number)
                elif type_output == 'cif':
                    name_file = 'structure_files/initial_structures/generated_structures/structure-' + "{:06d}".format(phase_number) + '.cif'
                struc.to_file(name_file, fmt=type_output)

                if text_output == True:
                    print(f'Crystal structure generated with the phase group {num_space_group + 1}!')

                phase_number = phase_number + 1
            except Exception:
                if text_output == True:
                    print(f'No compatibility of the stoichiometry with the phase group {num_space_group + 1}')
                else:
                    pass

    return phase_number


def relax_structure(relax_object, init_struc_path, relax_struc_path, max_steps, text_output, type_output, counter):
    """
    Relax an structure and get the energy

    Inputs:
        relax_object: the relaxer object by m3gnet module
        init_struc_path: the path of the structure to relax
        relax_struc_path: the path of the relaxed structure
        max_steps: number maximum of ionic steps
        text_output: to decide if we want output text in terminal or not
        type_output: the crystal type of output file, poscar or cif
        counter: a number counter for the phase
    """

    crystal_structure = Structure.from_file(init_struc_path)

    if text_output == True:
        relaxed_structure = relax_object.relax(crystal_structure, steps=max_steps, verbose=True)
    else:
        relaxed_structure = relax_object.relax(crystal_structure, steps=max_steps, verbose=False)

    final_structure = relaxed_structure['final_structure']

    final_structure.to(filename=relax_struc_path, fmt=type_output)

    final_energy = float(relaxed_structure['trajectory'].energies[-1]/len(crystal_structure))
    counter = counter + 1

    return final_energy, counter 


def eV_to_J(eV_energy):
    """
    This function returns energy in Joules from an energy in eV

    Inputs:
        eV_energy: energy in eV
    """

    return eV_energy*1.60218e-19


def J_to_eV(J_energy):
    """
    This function returns energy in eV from an energy in Joules

    Inputs:
        eV_energy: energy in eV
    """

    return J_energy*6.242e18


def A3_to_m3(A3_volume):
    """
    This function returns volume in m**3 from a volume in Angstrom**3

    Inputs:
        A3_volume: energy in Angstrom**3
    """

    return A3_volume*1e-30


def pressure_structure(relax_object, relax_struc_path, pressure_struc_path, pressure, number_volumes, min_alpha, max_alpha, text_output, type_output, counter):
    """
    Find the structure that minimizes enthalpy for a given pressure, and returnes the energy and alpha (where alpha is the proportional
    volume that that verifies the minimization)

    Inputs:
        relax_object: the relaxer object by m3gnet module
        relax_struc_path: the path of the relaxed structure
        pressure_struc_path: the path of the pressurized structure
        pressure: value of the pressure in Pa
        number_volumes: the number of different volumes we want compute in order to determine the enthalpy curve
        min_alpha: minimum proportional volume
        max_alpha: maximum proportional volume
        text_output: to decide if we want output text in terminal or not
        type_output: the crystal type of output file, poscar or cif
        counter: a number counter for the phase
    """

    relaxed_structure = Structure.from_file(relax_struc_path)

    volume = np.zeros(number_volumes)
    alpha_array = np.zeros(number_volumes)
    enthalpy = np.zeros(number_volumes)

    alpha = max_alpha

    for vol in range(number_volumes):
        distorted = deepcopy(relaxed_structure)
        distorted.scale_lattice(distorted.volume * alpha)

        volume[vol] = distorted.volume * alpha
        alpha_array[vol] = alpha

        if text_output == True:
            pressure_results = relax_object.relax(distorted, steps=1, verbose=True)
        else:
            pressure_results = relax_object.relax(distorted, steps=1, verbose=False)

        energy_volume = float(pressure_results['trajectory'].energies[0])

        enthalpy[vol] = J_to_eV(eV_to_J(energy_volume) + pressure*A3_to_m3(distorted.volume*alpha))

        alpha = alpha - (max_alpha - min_alpha)/number_volumes

    new_alpha_range = np.linspace(max_alpha, min_alpha, number_volumes*20)

    interpol = interp1d(alpha_array, enthalpy, kind='quadratic', fill_value='extrapolate')
    enthalpy_interpolated = interpol(new_alpha_range)

    alpha_min_enthalpy = new_alpha_range[np.argmin(enthalpy_interpolated)]

    pressurized_structure = deepcopy(relaxed_structure)
    pressurized_structure.scale_lattice(pressurized_structure.volume*alpha_min_enthalpy)

    pressurized_structure.to(filename=pressure_struc_path, fmt=type_output)

    if text_output == True:
        pressurized_energy = relax_object.relax(pressurized_structure, steps=1, verbose=True)
    else:
        pressurized_energy = relax_object.relax(pressurized_structure, steps=1, verbose=False)

    if alpha_min_enthalpy != new_alpha_range[-1]:
        final_energy = float(pressurized_energy['trajectory'].energies[0]/len(relaxed_structure))
    else:
        final_energy = 0

    counter = counter + 1

    return final_energy, counter, alpha_min_enthalpy


def write_in_file(file_object, structure_path, struc_num, energy, prec_spglib, type_output, counter):
    """
    Write the number, the file structure number, the energy per atom (in eV),
    and the phase group for a given structure in the energy file

    Inputs:
        file_object: the object of the file of energies
        structure_path: path to the structure of interest
        struc_num: the number of the structure
        energy: energy per atom (in eV) of the structure
        prec_spglib: precision for the phase determination with spglib
        type_output: the crystal type of output file, poscar or cif
        counter: a number counter for the phase
    """
    if type_output == 'poscar':
        structure = Poscar.from_file(structure_path).structure
    elif type_output == 'cif':
        parser = CifParser(structure_path)
        structure = parser.get_structures()[0]

    structure_symmetry = SpacegroupAnalyzer(structure=structure, symprec=prec_spglib)
    symmetry = structure_symmetry.get_symmetry_dataset()

    if type_output == 'poscar':
        file_object.write(f'{int(counter)}       POSCAR-{struc_num:06d}       {energy}       {symmetry["international"]} ({symmetry["number"]})\n')
    elif type_output == 'cif':
        file_object.write(f'{int(counter)}       structure-{struc_num:06d}.cif       {energy}       {symmetry["international"]} ({symmetry["number"]})\n')

    return


def write_in_file_pressure(file_object, structure_path, struc_num, energy, alpha, prec_spglib, type_output, counter):
    """
    Write the number, the file structure number, the energy per atom (in eV),
    and the phase group for a given structure in the energy file

    Inputs:
        file_object: the object of the file of energies
        structure_path: path to the structure of interest
        struc_num: the number of the structure
        energy: energy per atom (in eV) of the structure
        alpha: proportion of the volume of the pressurized phase
        prec_spglib: precision for the phase determination with spglib
        type_output: the crystal type of output file, poscar or cif
        counter: a number counter for the phase
    """
    if type_output == 'poscar':
        structure = Poscar.from_file(structure_path).structure
    elif type_output == 'cif':
        parser = CifParser(structure_path)
        structure = parser.get_structures()[0]

    structure_symmetry = SpacegroupAnalyzer(structure=structure, symprec=prec_spglib)
    symmetry = structure_symmetry.get_symmetry_dataset()

    if type_output == 'poscar':
        file_object.write(f'{int(counter)}       POSCAR-{struc_num:06d}       {energy}       {alpha}       {symmetry["international"]} ({symmetry["number"]})\n')
    elif type_output == 'cif':
        file_object.write(f'{int(counter)}       structure-{struc_num:06d}.cif       {energy}       {alpha}       {symmetry["international"]} ({symmetry["number"]})\n')

    return


def create_dir_generation(number_generation):
    """
    Creates the necessary directories to store all the files for each generation

    Inputs:
        number_generation: number of the actual generation
    """

    os.mkdir(f'structure_files/generation-{number_generation:03d}')
    os.mkdir(f'structure_files/generation-{number_generation:03d}/initial_structures')
    os.mkdir(f'structure_files/generation-{number_generation:03d}/distorted_structures')
    os.mkdir(f'structure_files/generation-{number_generation:03d}/relaxed_structures')
    os.mkdir(f'structure_files/generation-{number_generation:03d}/final_structures')

    return


def structure_distortion(file, max_displacement, type_output, final_path):
    """
    Distort a crystal structure

    Inputs:
        file: name or path of the file
        max_displacement: maximum displacement in a given direction (in Angstroms)
        type_output: the crystal type of output file, poscar or cif
        final_path: name and path of the final file
    """

    structure = Structure.from_file(file)

    displaced_structure = structure.copy()
    for site in displaced_structure:
        displacement_vector = np.random.uniform(-max_displacement, max_displacement, 3)
        site.coords = site.coords + displacement_vector

    displaced_structure.to(filename=final_path, fmt=type_output)

    return


def relax_structure_gen(relax_object, dist_struc_path, relax_struc_path, max_steps, text_output, type_output):
    """
    Relax an structure and get the energy for the generations part of the main code

    Inputs:
        relax_object: the relaxer object by m3gnet module
        dist_struc_path: the path of the distorted structure
        relax_struc_path: the path of the relaxed structure
        max_steps: number maximum of ionic steps
        text_output: to decide if we want output text in terminal or not
        type_output: the crystal type of output file, poscar or cif
    """

    crystal_structure = Structure.from_file(dist_struc_path)

    if text_output == True:
        relaxed_structure = relax_object.relax(crystal_structure, steps=max_steps, verbose=True)
    else:
        relaxed_structure = relax_object.relax(crystal_structure, steps=max_steps, verbose=False)

    final_structure = relaxed_structure['final_structure']

    final_structure.to(filename=relax_struc_path, fmt=type_output)

    final_energy = float(relaxed_structure['trajectory'].energies[-1]/len(crystal_structure))


    return final_energy


def write_in_file_gen(file_object, structure_path, struc_num, energy, prec_spglib, type_output, counter):
    """
    Write the number, the file structure number, the energy per atom (in eV),
    and the phase group for a given structure in the energy file

    Inputs:
        file_object: the object of the file of energies
        structure_path: path to the structure of interest
        struc_num: the number of the structure
        energy: energy per atom (in eV) of the structure
        prec_spglib: precision for the phase determination with spglib
        type_output: the crystal type of output file, poscar or cif
        counter: a number counter for the phase
    """

    if type_output == 'poscar':
        structure = Poscar.from_file(structure_path).structure
    elif type_output == 'cif':
        parser = CifParser(structure_path)
        structure = parser.get_structures()[0]

    structure_symmetry = SpacegroupAnalyzer(structure=structure, symprec=prec_spglib)
    symmetry = structure_symmetry.get_symmetry_dataset()

    file_object.write(f'{int(counter)}       {struc_num}       {energy}       {symmetry["international"]} ({symmetry["number"]})\n')


    return
