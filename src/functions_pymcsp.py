# Pol Benítez Colominas, Universitat Politècnica de Catalunya
# September 2023 - January 2025
# Version 1.0

# Functions file

import os
import sys
import shutil
import time
from copy import deepcopy
import csv
import ast

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy import signal
from scipy import integrate

from pymatgen.core import Structure
from pymatgen.io.vasp import Poscar
from pymatgen.io.cif import CifParser
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.ase import AseAtomsAdaptor

from pyxtal import pyxtal

import matgl
from matgl.ext.ase import Relaxer, PESCalculator

from terminal_outputs import *

import warnings

warnings.simplefilter("ignore") # To suppress warnings for clearer output


def read_variables_csp(path_file):
    """
    Read the variables from the inputs_csp file and returns an array with their values

    Inputs:
        path_file: path of the inputs file
    """

    inputs = open(path_file, "r")

    variables = [None]*22

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

    return variables


def read_variables_diffraction(path_file):
    """
    Read the variables from the inputs_diffraction file and returns an array with their values

    Inputs:
        path_file: path of the inputs file
    """

    inputs = open(path_file, "r")

    variables = [None]*21

    for it in range(12):
        inputs.readline()

    for it in range(len(variables)):
        line = inputs.readline()

        if line.split()[2] == 'True':
            variables[it] = True
        elif line.split()[2] == 'False':
            variables[it] = False
        else:
            variables[it] = line.split()[2]

    inputs.close()

    return variables


def generate_phases(stoi_dict, phase_number, dim, restricted, restricted_list, atoms_arr, text_output, type_output):
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

        if restricted == False:
            for num_space_group in range(num_groups[dim]):
                try:
                    struc.from_random(dim, num_space_group+1, atoms_arr, value)
                    if type_output == 'poscar':
                        name_file = 'structure_files/generated_structures/POSCAR-' + "{:06d}".format(phase_number)
                    elif type_output == 'cif':
                        name_file = 'structure_files/generated_structures/structure-' + "{:06d}".format(phase_number) + '.cif'
                    struc.to_file(name_file, fmt=type_output)

                    if text_output == True:
                        print(f'Crystal structure generated with the phase group {num_space_group + 1}!')

                    phase_number = phase_number + 1
                except Exception:
                    if text_output == True:
                        print(f'No compatibility of the stoichiometry with the phase group {num_space_group + 1}')
                    else:
                        pass
        elif restricted == True:
            for num_space_group in restricted_list:
                try:
                    struc.from_random(dim, num_space_group, atoms_arr, value)
                    if type_output == 'poscar':
                        name_file = 'structure_files/generated_structures/POSCAR-' + "{:06d}".format(phase_number)
                    elif type_output == 'cif':
                        name_file = 'structure_files/generated_structures/structure-' + "{:06d}".format(phase_number) + '.cif'
                    struc.to_file(name_file, fmt=type_output)

                    if text_output == True:
                        print(f'Crystal structure generated with the phase group {num_space_group}!')

                    phase_number = phase_number + 1
                except Exception:
                    if text_output == True:
                        print(f'No compatibility of the stoichiometry with the phase group {num_space_group}')
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


def pressure_structure(calc_object, ase_adaptor_object, relax_struc_path, pressure_struc_path, pressure, number_volumes, min_alpha, max_alpha, text_output, type_output, counter):
    """
    Find the structure that minimizes enthalpy for a given pressure, and returnes the energy and alpha (where alpha is the proportional
    volume that that verifies the minimization)

    Inputs:
        calc_object: the single point energy calculation object by m3gnet
        ase_adaptor_object: ase adaptator object necessary for single point energy calculation
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

    num_atoms = relaxed_structure.num_sites

    volume = np.zeros(number_volumes)
    alpha_array = np.zeros(number_volumes)
    enthalpy = np.zeros(number_volumes)

    alpha = max_alpha

    for vol in range(number_volumes):
        distorted = deepcopy(relaxed_structure)
        distorted.scale_lattice(distorted.volume * alpha)

        volume[vol] = distorted.volume * alpha
        alpha_array[vol] = alpha

        atoms = ase_adaptor_object.get_atoms(distorted)
        atoms.set_calculator(calc_object)

        energy_volume = float(atoms.get_potential_energy())

        enthalpy[vol] = J_to_eV(eV_to_J(energy_volume) + pressure*A3_to_m3(distorted.volume*alpha))

        alpha = alpha - (max_alpha - min_alpha)/number_volumes

    new_alpha_range = np.linspace(max_alpha, min_alpha, number_volumes*20)

    interpol = interp1d(alpha_array, enthalpy, kind='quadratic', fill_value='extrapolate')
    enthalpy_interpolated = interpol(new_alpha_range)

    alpha_min_enthalpy = new_alpha_range[np.argmin(enthalpy_interpolated)]

    pressurized_structure = deepcopy(relaxed_structure)
    pressurized_structure.scale_lattice(pressurized_structure.volume*alpha_min_enthalpy)

    pressurized_structure.to(filename=pressure_struc_path, fmt=type_output)

    atoms = ase_adaptor_object.get_atoms(pressurized_structure)
    atoms.set_calculator(calc_object)

    if alpha_min_enthalpy != new_alpha_range[-1]:
        final_energy = float(atoms.get_potential_energy())
        final_enthalpy = final_energy + J_to_eV((pressure*alpha_min_enthalpy*relaxed_structure.volume*(1e-30))/num_atoms)
    else:
        final_energy = 0
        final_enthalpy = 0

    counter = counter + 1

    return final_enthalpy, final_energy, counter, alpha_min_enthalpy


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


def write_in_file_pressure(file_object, structure_path, struc_num, enthalpy, energy, alpha, prec_spglib, type_output, counter):
    """
    Write the number, the file structure number, the energy per atom (in eV),
    and the phase group for a given structure in the energy file

    Inputs:
        file_object: the object of the file of energies
        structure_path: path to the structure of interest
        struc_num: the number of the structure
        enthalpy: determined enthalpy per atom
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
        file_object.write(f'{int(counter)}       POSCAR-{struc_num:06d}       {enthalpy}       {energy}       {alpha}       {symmetry["international"]} ({symmetry["number"]})\n')
    elif type_output == 'cif':
        file_object.write(f'{int(counter)}       structure-{struc_num:06d}.cif       {enthalpy}       {energy}       {alpha}       {symmetry["international"]} ({symmetry["number"]})\n')

    return


def deduplicate_file(path_energy_file, struc_path, num_structures, type_output, prec_group):
    """
    Creates a file with deduplicated energy phases

    Inputs:
        path_energy_file: path to the energy ranking of structures file
        struc_path: path to the directory with the structures
        num_structures: total number of structures (in the directory and the energy ranking file (same))
        type_output: the crystal type of output file, poscar or cif
        prec_group: precision for the phase determination with spglib
    """

    energy_threshold = 0.0005 # 0.5 meV/atom, typical computational resolution of DFT

    non_considered_structures = list(range(num_structures))

    energy_file = open(path_energy_file, 'r')
    energy_file_deduplicate = open(struc_path + 'energy_deduplicate.txt', 'w')

    line_file = energy_file.readline()
    energy_file_deduplicate.write(line_file)

    num_new_struc = 0
    for struc in range(num_structures):
        line_file = energy_file.readline()

        file_name = line_file.split()[1]
        energy_struc = float(line_file.split()[2])
        space_group = line_file.split()[3] + ' ' + line_file.split()[4]

        if (struc in non_considered_structures) and (struc != (num_structures - 1)):
            non_considered_structures.remove(struc)
            energy_file2 = open(path_energy_file, 'r')
            for iter in range(struc + 2):
                energy_file2.readline()
            
            finish_condition = False
            while_it = 0
            while finish_condition == False:
                line_file2 = energy_file2.readline()

                file_name2 = line_file2.split()[1]
                energy_struc2 = float(line_file2.split()[2])
                space_group2 = line_file2.split()[3] + ' ' + line_file2.split()[4]

                while_it = while_it + 1

                if abs(energy_struc - energy_struc2) <= energy_threshold:
                    if energy_struc >= energy_struc2:
                        if type_output == 'poscar':
                            structure = Poscar.from_file(struc_path + '/' + file_name).structure
                            structure2 = Poscar.from_file(struc_path + '/' + file_name2).structure
                        elif type_output == 'cif':
                            parser = CifParser(struc_path + '/' + file_name)
                            structure = parser.get_structures()[0]
                            parser2 = CifParser(struc_path + '/' + file_name2)
                            structure2 = parser2.get_structures()[0]

                        num_atoms = structure.num_sites
                        num_atoms2 = structure2.num_sites

                        if num_atoms > num_atoms2:
                            file_name = file_name2
                            energy_struc = energy_struc2
                            space_group = space_group2
                        elif num_atoms == num_atoms2:
                            structure_symmetry = SpacegroupAnalyzer(structure=structure, symprec=prec_group)
                            symmetry = structure_symmetry.get_space_group_number()
                            structure_symmetry2 = SpacegroupAnalyzer(structure=structure2, symprec=prec_group)
                            symmetry2 = structure_symmetry2.get_space_group_number()

                            if symmetry < symmetry2:
                                file_name = file_name2
                                energy_struc = energy_struc2
                                space_group = space_group2
                                
                    non_considered_structures.remove(struc + while_it)
                else:
                    finish_condition = True
                    non_considered_structures.remove(struc + while_it)

                if (struc + while_it) == num_structures - 1:
                    finish_condition = True
            
            num_new_struc = num_new_struc + 1

            energy_file_deduplicate.write(f'{num_new_struc}       {file_name}       {energy_struc}       {space_group}\n')

        elif (struc in non_considered_structures) and (struc == (num_structures - 1)):
            num_new_struc = num_new_struc + 1

            energy_file_deduplicate.write(f'{num_new_struc}       {file_name}       {energy_struc}       {space_group}\n')

    energy_file.close()
    energy_file_deduplicate.close()

    return


def deduplicate_file_pressure(path_energy_file, struc_path, num_structures, type_output, prec_group):
    """
    Creates a file with deduplicated energy phases for pressurized results

    Inputs:
        path_energy_file: path to the energy ranking of structures file
        struc_path: path to the directory with the structures
        num_structures: total number of structures (in the directory and the energy ranking file (same))
        type_output: the crystal type of output file, poscar or cif
        prec_group: precision for the phase determination with spglib
    """

    energy_threshold = 0.0005 # 0.5 meV/atom, typical computational resolution of DFT

    non_considered_structures = list(range(num_structures))

    energy_file = open(path_energy_file, 'r')
    energy_file_deduplicate = open(struc_path + 'energy_deduplicate.txt', 'w')

    line_file = energy_file.readline()
    energy_file_deduplicate.write(line_file)

    num_new_struc = 0
    for struc in range(num_structures):
        line_file = energy_file.readline()

        file_name = line_file.split()[1]
        enthalpy_struc = float(line_file.split()[2])
        energy_struc = float(line_file.split()[3])
        alpha_coef = float(line_file.split()[4])
        space_group = line_file.split()[5] + ' ' + line_file.split()[6]

        if (struc in non_considered_structures) and (struc != (num_structures - 1)):
            non_considered_structures.remove(struc)
            energy_file2 = open(path_energy_file, 'r')
            for iter in range(struc + 2):
                energy_file2.readline()
            
            finish_condition = False
            while_it = 0
            while finish_condition == False:
                line_file2 = energy_file2.readline()

                file_name2 = line_file2.split()[1]
                enthalpy_struc2 = float(line_file2.split()[2])
                energy_struc2 = float(line_file2.split()[3])
                alpha_coef2 = float(line_file2.split()[4])
                space_group2 = line_file2.split()[5] + ' ' + line_file2.split()[6]

                while_it = while_it + 1

                if abs(energy_struc - energy_struc2) <= energy_threshold:
                    if energy_struc >= energy_struc2:
                        if type_output == 'poscar':
                            structure = Poscar.from_file(struc_path + '/' + file_name).structure
                            structure2 = Poscar.from_file(struc_path + '/' + file_name2).structure
                        elif type_output == 'cif':
                            parser = CifParser(struc_path + '/' + file_name)
                            structure = parser.get_structures()[0]
                            parser2 = CifParser(struc_path + '/' + file_name2)
                            structure2 = parser2.get_structures()[0]

                        num_atoms = structure.num_sites
                        num_atoms2 = structure2.num_sites

                        if num_atoms > num_atoms2:
                            file_name = file_name2
                            enthalpy_struc = enthalpy_struc2
                            energy_struc = energy_struc2
                            alpha_coef = alpha_coef2
                            space_group = space_group2
                        elif num_atoms == num_atoms2:
                            structure_symmetry = SpacegroupAnalyzer(structure=structure, symprec=prec_group)
                            symmetry = structure_symmetry.get_space_group_number()
                            structure_symmetry2 = SpacegroupAnalyzer(structure=structure2, symprec=prec_group)
                            symmetry2 = structure_symmetry2.get_space_group_number()

                            if symmetry < symmetry2:
                                file_name = file_name2
                                enthalpy_struc = enthalpy_struc2
                                energy_struc = energy_struc2
                                alpha_coef = alpha_coef2
                                space_group = space_group2
                                
                    non_considered_structures.remove(struc + while_it)
                else:
                    finish_condition = True
                    non_considered_structures.remove(struc + while_it)

                if (struc + while_it) == num_structures - 1:
                    finish_condition = True
            
            num_new_struc = num_new_struc + 1

            energy_file_deduplicate.write(f'{num_new_struc}       {file_name}       {enthalpy_struc}       {energy_struc}       {alpha_coef}       {space_group}\n')

        elif (struc in non_considered_structures) and (struc == (num_structures - 1)):
            num_new_struc = num_new_struc + 1

            energy_file_deduplicate.write(f'{num_new_struc}       {file_name}       {enthalpy_struc}       {energy_struc}       {alpha_coef}       {space_group}\n')

    energy_file.close()
    energy_file_deduplicate.close()

    return


def create_dir_generation(number_generation):
    """
    Creates the necessary directories to store all the files for each generation

    Inputs:
        number_generation: number of the actual generation
    """

    os.mkdir(f'structure_files/generations/generation-{number_generation:03d}')
    os.mkdir(f'structure_files/generations/generation-{number_generation:03d}/initial_structures')
    os.mkdir(f'structure_files/generations/generation-{number_generation:03d}/distorted_structures')
    os.mkdir(f'structure_files/generations/generation-{number_generation:03d}/relaxed_structures')
    os.mkdir(f'structure_files/generations/generation-{number_generation:03d}/final_structures')

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


def curve_from_peaks(prev_peaks_x, prev_peaks_y, min_x, max_x, num_points, kind_interpolation, peak_width):
    """
    Returns a curve from a given set of peaks (neglects those peaks that are small, smaller that 1e-1)

    Inputs:
        peaks_x: value of the peaks in x range
        peaks_y: value of the peaks in y range
        min_x: minimum value of the x range of the interpolation
        max_x: maximum value of the x range of the interpolation
        num_points: number of points of the interpolation
        kind_interpolation: type of interpolation performed
        peak_width: width of the peak
    """
    peaks_x = []
    peaks_y = []
    for peak in range(len(prev_peaks_y)):
        if prev_peaks_y[peak] >= 1e-1:
            peaks_x.append(prev_peaks_x[peak])
            peaks_y.append(prev_peaks_y[peak])

    curve_x = np.zeros(3*len(peaks_x) + 2)
    curve_y = np.zeros(3*len(peaks_x) + 2)

    counter_peaks = 0
    for iteration in range(len(curve_x)):
        if (iteration != 0) and (iteration != (len(curve_x) - 1)):
            if (iteration + 1) % 3 == 0:
                curve_x[iteration] = peaks_x[counter_peaks]
                curve_y[iteration] = peaks_y[counter_peaks]

                counter_peaks = counter_peaks + 1
            else:
                if (iteration) % 3 == 0:
                    curve_x[iteration] = peaks_x[counter_peaks - 1] + peak_width
                elif (iteration + 2) % 3 == 0:
                    curve_x[iteration] = peaks_x[counter_peaks] - peak_width

                curve_y[iteration] = 0 

    curve_x[0] = curve_x[1] - peak_width
    curve_y[0] = 0
    curve_x[-1] = curve_x[-2] + peak_width
    curve_y[-1] = 0

    range_x = np.linspace(min_x, max_x, num_points)
    interpol = interp1d(curve_x, curve_y, kind=kind_interpolation, fill_value='extrapolate')
    range_y = interpol(range_x)

    return range_x, range_y


def normalize_curve(curve_x, curve_y):
    """
    Retunrs a normalized curve

    Inputs:
        curve_x: range of values in x
        curve_y: range of values in y
    """

    normalization_ct = 1 / (integrate.trapezoid(curve_y, curve_x))

    range_y = np.zeros(len(curve_y))
    range_y[:] = normalization_ct*curve_y[:]

    return range_y


def peaks_experimental(curve_x, curve_y, prominance, width):
    """
    Find the peaks of experimental diffractogram curve

    Inputs:
        curve_x: range of values in x
        curve_y: range of values in y
        prominance, width: https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.find_peaks.html
    """

    peaks, _ = signal.find_peaks(curve_y, prominence=prominance, width=width)

    peak_x = np.zeros(len(peaks))
    peak_y = np.zeros(len(peaks))
    for x in range(len(peaks)):
        peak_x[x] = curve_x[peaks[x]]
        peak_y[x] = curve_y[peaks[x]]

    return peak_x, peak_y


def compute_loss_factor(range_x, intensity_theoretical, intesity_experimental, num_peaks_theo, num_peaks_exp, coef_quad):
    """
    Computes the loss factor of similarity between theoretical and experimental curves
    bigger values represent less similarity. This function also penalizes a big difference
    in the number of points

    Inputs:
        range_x: range of values in x
        intensity_theoretical: intesity values of the theoretical phase
        intensity_experimental: intesity values of the obtained experimental diffractogram
        num_peaks_theo: number of peaks in the theoretical diffraction curve
        num_peaks_exp: number of peaks in the experimental diffraction curve
        coef_quad: quadratic coefficient for difference of peak number
    """

    function_to_integrate = np.zeros(len(range_x))
    function_to_integrate[:] = np.abs(intesity_experimental[:] - intensity_theoretical[:])**2

    T_result = integrate.trapezoid(function_to_integrate, range_x) + coef_quad*np.abs(num_peaks_theo - num_peaks_exp)**2

    return T_result


def compute_loss_factor_minim(range_x, intensity_theoretical, intesity_experimental):
    """
    Computes the loss factor of similarity between theoretical and experimental curves
    bigger values represent less similarity, without taking into account the difference 
    in the number of peaks. This is used in minimize_loss_vol_scale function

    Inputs:
        range_x: range of values in x
        intensity_theoretical: intesity values of the theoretical phase
        intensity_experimental: intesity values of the obtained experimental diffractogram
    """

    function_to_integrate = np.zeros(len(range_x))
    function_to_integrate[:] = np.abs(intesity_experimental[:] - intensity_theoretical[:])**2

    T_result = integrate.trapezoid(function_to_integrate, range_x)

    return T_result


def minimize_loss_vol_scale(structure, calculator, min_x, max_x, num_points, kind_interpolation, peak_width, prop_vol, num_vols, intesity_experimental):
    """
    Scales the volume of the structure in the desired interval and finds the structure with minimum loss factor

    Inputs:
        structure: structure pymatgen object to study
        calculator: calculator object
        min_x: minimum value of the x range of the interpolation
        max_x: maximum value of the x range of the interpolation
        num_points: number of points of the interpolation
        kind_interpolation: type of interpolation performed
        peak_width: width of the peak
        prop_vol: proportion the volume will be increase and decrease, example: if prop_vol=0.05 it will be considered volumes between 0.95 an 1.05 the original volume
        num_vols: number of volumes to study between (1 - prop_vol) and (1 + prop_vol)
        intensity_experimental: intesity values of the obtained experimental diffractogram
    """

    vol_array = np.linspace(1 - prop_vol, 1 + prop_vol, num_vols)
    loss_array = []

    for vol in vol_array:
        distorted = deepcopy(structure)
        distorted.scale_lattice(distorted.volume * vol)

        pattern = calculator.get_pattern(distorted)
        two_theta_peaks = pattern.x
        intensity_peaks = pattern.y

        twotheta, intensity = curve_from_peaks(two_theta_peaks, intensity_peaks, min_x, max_x, num_points, kind_interpolation, peak_width)

        normalized_intensity = normalize_curve(twotheta, intensity)

        loss_factor = compute_loss_factor_minim(twotheta, normalized_intensity, intesity_experimental)
        loss_array.append(loss_factor)

    min_loss_index = loss_array.index(min(loss_array))
    vol_to_minimize = vol_array[min_loss_index]

    return vol_to_minimize

def remove_substrate(exp_peaks_x, exp_peaks_y, subs_peaks_x, tolerance):
    """
    Returns the experimental peaks without the peaks related with the substrate, to do this
    the distance between two peaks is computed and if the value it is below a given tolerance
    that peak will be discarted

    Inputs:
        exp_peaks_x: 2theta value of the experimental peaks
        exp_peaks_y: intensity value of the experimental peaks
        subs_peaks_x: 2theta value of the theoretical substrate peaks
        tolerance: threshold value to discard a peak
    """

    clean_peaks_x = []
    clean_peaks_y = []

    for iter in range(len(exp_peaks_x)):
        substrate_peak = False
        for iter2 in range(len(subs_peaks_x)):
            peak_distance = np.abs(exp_peaks_x[iter] - subs_peaks_x[iter2])
            if peak_distance <= tolerance:
                substrate_peak = True
        
        if substrate_peak == False:
            clean_peaks_x.append(exp_peaks_x[iter])
            clean_peaks_y.append(exp_peaks_y[iter])

    return clean_peaks_x, clean_peaks_y


def plot_exp_tuned(prominance, width, exp_path):
    """"
    Find the peaks of the experimental curve with the provided parameters and show them in a plot

    Inputs:
        prominance: how much peak stands out respect the baseline
        width: width of the peaks
        inputs_path: path to inputs_diffraction file
    """

    exp_2theta = []
    exp_int = []

    with open(exp_path, 'r') as file:
        reader = csv.reader(file)

        has_header = csv.Sniffer().has_header(file.read(1024))
        file.seek(0)
        
        if has_header:
            next(reader) 

        for row in reader:
            exp_2theta.append(float(row[0]))  
            exp_int.append(float(row[1]))  

    exp_2theta = np.array(exp_2theta)
    exp_int = np.array(exp_int)

    exp_peaks_x, exp_peaks_y = peaks_experimental(exp_2theta, exp_int, prominance, width)

    plt.figure()
    plt.plot(exp_2theta, exp_int, label='Experimental curve')
    plt.plot(exp_peaks_x, exp_peaks_y, marker='o', linestyle='', label='Found peaks')
    plt.legend()
    plt.show()

    return


def relaxer_function(retrained, retrained_path):
    """
    Calls M3GNet relaxer function

    Inputs:
        retreined: boolean, True if the user provides a retrained model path
        retreined_path: path to the retreined model for M3GNet
    """

    if retrained == False:
        pot = matgl.load_model("M3GNet-MP-2021.2.8-PES")
    else:
        pot = matgl.load_model(retrained_path)

    relaxer = Relaxer(potential=pot)

    return relaxer


def energy_function(retrained, retrained_path):
    """
    Calls M3GNet single point calculation function

    Inputs:
        retreined: boolean, True if the user provides a retrained model path
        retreined_path: path to the retreined model for M3GNet
    """

    if retrained == False:
        pot = matgl.load_model("M3GNet-MP-2021.2.8-PES")
    else:
        pot = matgl.load_model(retrained_path)

    calc = PESCalculator(potential=pot)

    return calc


def csp_study(inputs):
    """
    Performs crystal structure prediction

    Inputs:
        inputs: array with all the input variables from inputs_csp file
    """

    start_time = time.time()

    dimension = int(inputs[0])
    atoms = inputs[1]
    stoichiometry = inputs[2]
    num_max_atoms = int(inputs[3])
    num_same_phase = int(inputs[4])
    max_ionic_steps = int(inputs[5])
    print_terminal_outputs = inputs[10]
    save_log_file = inputs[11]
    structure_file = inputs[12]
    prec_group_det = float(inputs[13])
    retrain = inputs[18]
    retrain_path = inputs[19]
    restricted = inputs[20]
    restricted_list = ast.literal_eval(inputs[21])


    ### Save the log file
        
    log_file = 'output_csp.log'

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

        simulation_information(atoms, stoichiometry, num_max_atoms)


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

    os.mkdir('structure_files')
    os.mkdir('structure_files/generated_structures/')

    num_phase = 1
    for num_phase_group in range(num_same_phase):
        total_num_phase = generate_phases(stoichiometry_dict, num_phase, dimension, restricted,
                                        restricted_list, atoms, print_terminal_outputs, structure_file)
        
        num_phase = total_num_phase

    if print_terminal_outputs == True:
        struc_gen_end()

    
    #### Relax the first structures
    relaxer = relaxer_function(retrained=retrain, retrained_path=retrain_path)

    dir_path = 'structure_files/generated_structures/'

    num_structures = 0

    for path in os.listdir(dir_path):
        if os.path.isfile(os.path.join(dir_path, path)):
            num_structures = num_structures + 1

    if print_terminal_outputs == True:
        print(f'Number of structures to relax: {num_structures}')

    os.mkdir('structure_files/relaxed_structures/')

    phase_energy_array = np.zeros((num_structures,2))
    count = 0

    if print_terminal_outputs == True:
        relax_ini()

    for num_struc in range(num_structures):
        if structure_file == 'poscar':
            initial_path = 'structure_files/generated_structures/POSCAR-' + "{:06d}".format(num_struc + 1)
            relaxed_path = 'structure_files/relaxed_structures/POSCAR-' + "{:06d}".format(num_struc + 1)
        elif structure_file == 'cif':
            initial_path = 'structure_files/generated_structures/structure-' + "{:06d}".format(num_struc + 1) + '.cif'
            relaxed_path = 'structure_files/relaxed_structures/structure-' + "{:06d}".format(num_struc + 1) + '.cif'

        if print_terminal_outputs == True:
            if structure_file == 'poscar':
                struc_name = 'POSCAR-' + "{:06d}".format(num_struc + 1)
            elif structure_file == 'cif':
                struc_name = 'structure-' + "{:06d}".format(num_struc + 1) + '.cif'
            relax_struc(struc_name)

        relax_energy, count = relax_structure(relaxer, initial_path, relaxed_path, max_ionic_steps, 
                                            print_terminal_outputs, structure_file, count)
        
        phase_energy_array[count - 1, 0] = int(count)
        phase_energy_array[count - 1, 1] = relax_energy

    if print_terminal_outputs == True:
        relax_end()

    sorted_indices = np.argsort(phase_energy_array[:,1])
    phase_energy_sorted = phase_energy_array[sorted_indices]

    file_energy = open('structure_files/relaxed_structures/energy_ranking.txt', 'w')
    if structure_file == 'poscar':
        file_energy.write('#       POSCAR-num       energy per atom (eV)       Space Group (Hermann-Mauguin)\n')
    elif structure_file == 'cif':
        file_energy.write('#       structure-num.cif       energy per atom (eV)       Space Group (Hermann-Mauguin)\n')

    count = 1
    for num_struc in range(num_structures):
        struc_number = int(phase_energy_sorted[num_struc,0])
        struc_energy = phase_energy_sorted[num_struc,1]
        if structure_file == 'poscar':
            phase_file = 'structure_files/relaxed_structures/POSCAR-' + "{:06d}".format(struc_number)
        elif structure_file == 'cif':
            phase_file = 'structure_files/relaxed_structures/structure-' + "{:06d}".format(struc_number) + '.cif'

        write_in_file(file_energy, phase_file, struc_number, struc_energy, 
                    prec_group_det, structure_file, count)

        count = count + 1

    file_energy.close()

    deduplicate_file('structure_files/relaxed_structures/energy_ranking.txt', 'structure_files/relaxed_structures/', 
                     num_structures, structure_file, prec_group_det)

    end_time = time.time()
    elapsed_time = end_time - start_time

    if print_terminal_outputs == True:
        print(f'Elapsed time: {elapsed_time:.3f} seconds')


    ### Close the log file
        
    if save_log_file == True:
        log_file_handle.close()
        sys.stdout = sys.__stdout__

    return


def pressure_computations(inputs):
    """
    Performs pressure computations for an already found structures

    Inputs:
        inputs: array with all the input variables from inputs_csp file
    """

    start_time = time.time()

    pressure = float(inputs[6])
    num_volumes = int(inputs[7])
    minimum_volume = float(inputs[8]) 
    maximum_volume = float(inputs[9])
    print_terminal_outputs = inputs[10]
    save_log_file = inputs[11]
    structure_file = inputs[12]
    prec_group_det = float(inputs[13])
    struc_path = inputs[17]
    retrain = inputs[18]
    retrain_path = inputs[19]

    ### Save the log file
        
    log_file = 'output_pressure.log'

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

    ### Start the pressure computation

    if print_terminal_outputs == True:
        initial_message()

        simulation_information_pressure(pressure)

    num_structures = 0

    for path in os.listdir(struc_path):
        if os.path.isfile(os.path.join(struc_path, path)):
            num_structures = num_structures + 1

    if print_terminal_outputs == True:
        print(f'Number of structures to study with pressure conditions: {num_structures - 2}')

    if os.path.exists('structure_files/pressure_structures/'):
        shutil.rmtree('structure_files/pressure_structures/')
    os.mkdir('structure_files/pressure_structures/')

    pressure_energy_array = np.zeros((num_structures - 2, 4))
    count = 0

    if print_terminal_outputs == True:
        press_ini()

    calc = energy_function(retrained=retrain, retrained_path=retrain_path)
    ase_adaptor = AseAtomsAdaptor()

    for num_struc in range(num_structures - 2):
        if structure_file == 'poscar':
            relaxed_path = struc_path + '/POSCAR-' + "{:06d}".format(num_struc + 1)
            pressure_path = 'structure_files/pressure_structures/POSCAR-' + "{:06d}".format(num_struc + 1)
        elif structure_file == 'cif':
            relaxed_path = struc_path +  '/structure-' + "{:06d}".format(num_struc + 1) + '.cif'
            pressure_path = 'structure_files/pressure_structures/structure-' + "{:06d}".format(num_struc + 1) + '.cif'

        if print_terminal_outputs == True:
            if structure_file == 'poscar':
                struc_name = 'POSCAR-' + "{:06d}".format(num_struc + 1)
            elif structure_file == 'cif':
                struc_name = 'structure-' + "{:06d}".format(num_struc + 1) + '.cif'
            press_struc(struc_name)

        pressure_enthalpy, pressure_energy, count, alpha = pressure_structure(calc, ase_adaptor, relaxed_path, pressure_path, pressure, num_volumes, minimum_volume,
                                                           maximum_volume, print_terminal_outputs, structure_file, count)
            
        pressure_energy_array[count - 1, 0] = int(count)
        pressure_energy_array[count - 1, 1] = pressure_energy
        pressure_energy_array[count - 1, 2] = alpha
        pressure_energy_array[count - 1, 3] = pressure_enthalpy

    if print_terminal_outputs == True:
        press_end()

    sorted_indices_pressure = np.argsort(pressure_energy_array[:,3])
    phase_energy_sorted_pressure = pressure_energy_array[sorted_indices_pressure]

    file_energy = open('structure_files/pressure_structures/energy_ranking_pressure.txt', 'w')
    if structure_file == 'poscar':
        file_energy.write('#       POSCAR-num       enthalpy per atom (eV)       energy per atom (eV)       alpha       Space Group (Hermann-Mauguin)\n')
    elif structure_file == 'cif':
        file_energy.write('#       structure-num.cif       enthalpy per atom (eV)       energy per atom (eV)       alpha       Space Group (Hermann-Mauguin)\n')

    count = 1
    for num_struc in range(num_structures - 2):
        struc_number = int(phase_energy_sorted_pressure[num_struc,0])
        struc_energy = phase_energy_sorted_pressure[num_struc,1]
        struc_alpha = phase_energy_sorted_pressure[num_struc,2]
        struc_enthalpy = phase_energy_sorted_pressure[num_struc,3]
        if structure_file == 'poscar':
            phase_file = 'structure_files/pressure_structures/POSCAR-' + "{:06d}".format(struc_number)
        elif structure_file == 'cif':
            phase_file = 'structure_files/pressure_structures/structure-' + "{:06d}".format(struc_number) + '.cif'

        if struc_energy == 0:
            struc_energy = 'Error'
            struc_enthalpy = 'Error'
            struc_alpha = 'Not minimized'

        write_in_file_pressure(file_energy, phase_file, struc_number, struc_enthalpy, struc_energy, struc_alpha,
                        prec_group_det, structure_file, count)

        count = count + 1

    file_energy.close()

    deduplicate_file_pressure('structure_files/pressure_structures/energy_ranking_pressure.txt', 'structure_files/pressure_structures/', 
                     num_structures - 2, structure_file, prec_group_det)

    end_time = time.time()
    elapsed_time = end_time - start_time

    if print_terminal_outputs == True:
        print(f'Elapsed time: {elapsed_time:.3f} seconds')


    ### Close the log file
        
    if save_log_file == True:
        log_file_handle.close()
        sys.stdout = sys.__stdout__

    return


def generations_loop(inputs):
    """
    Performs generations loop for an already found structures

    Inputs:
        inputs: array with all the input variables from inputs_csp file
    """

    start_time = time.time()

    max_ionic_steps = int(inputs[5])
    print_terminal_outputs = inputs[10]
    save_log_file = inputs[11]
    structure_file = inputs[12]
    prec_group_det = float(inputs[13])
    num_generations = int(inputs[14])
    surviving_phases = float(inputs[15])
    max_disp = float(inputs[16])
    struc_path = inputs[17]
    retrain = inputs[18]
    retrain_path = inputs[19]

    ### Save the log file
        
    log_file = 'output_generations.log'

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

    ### Start the generation loop

    files_in_dir = os.listdir(struc_path)

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

    if os.path.exists('structure_files/generations/'):
        shutil.rmtree('structure_files/generations/')
    os.mkdir('structure_files/generations/')

    relaxer = relaxer_function(retrained=retrain, retrained_path=retrain_path)

    if print_terminal_outputs == True:
        gen_ini()

    for num_gen in range(num_generations):

        if print_terminal_outputs == True:
            gen_ini_actual(num_gen + 1)

        create_dir_generation(num_gen + 1)

        name_dir_initial_struc = {
            'initial': struc_path,
            'not_initial': f'structure_files/generations/generation-{num_gen:03d}/final_structures'
        }

        path_destination = f'structure_files/generations/generation-{num_gen + 1:03d}/initial_structures/'
        if structure_file == 'poscar':
            file_prefix = 'POSCAR-'
        elif structure_file == 'cif':
            file_prefix = 'structure-'

        if  num_gen == 0:
            path_initial_struc = name_dir_initial_struc['initial']

            name_selected_structures = []

            energy_file = open(struc_path + '/energy_ranking.txt', "r")

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
            path_struc_file = f'structure_files/generations/generation-{num_gen + 1:03d}/initial_structures/' + struc_file
            path_dist_file = f'structure_files/generations/generation-{num_gen + 1:03d}/distorted_structures/' + struc_file

            structure_distortion(path_struc_file, max_disp, structure_file, path_dist_file)

        file_energy = open(f'structure_files/generations/generation-{num_gen + 1:03d}/relaxed_structures/energy_ranking.txt', 'w')
        if structure_file == 'poscar':
            file_energy.write('#       POSCAR-num       energy per atom (eV)       Space Group (Hermann-Mauguin)\n')
        elif structure_file == 'cif':
            file_energy.write('#       structure-num.cif       energy per atom (eV)       Space Group (Hermann-Mauguin)\n')
        count = 1

        if print_terminal_outputs == True:
            relax_ini()

        for struc_file in name_selected_structures:
            path_dist_file = f'structure_files/generations/generation-{num_gen + 1:03d}/distorted_structures/' + struc_file
            path_relax_file = f'structure_files/generations/generation-{num_gen + 1:03d}/relaxed_structures/' + struc_file

            if print_terminal_outputs == True:
                relax_struc(struc_file)

            relaxed_energy = relax_structure_gen(relaxer, path_dist_file, path_relax_file, max_ionic_steps,
                                                    print_terminal_outputs, structure_file)
                
            write_in_file_gen(file_energy, path_relax_file, struc_file, relaxed_energy, 
                                prec_group_det, structure_file, count)
                
            count = count + 1

        if print_terminal_outputs == True:
            relax_end()

        file_energy.close()

        path_previous_energy_file = {
            'initial': struc_path + '/energy_ranking.txt',
            'not_initial': f'structure_files/generations/generation-{num_gen:03d}/final_structures/energy_ranking.txt'
        }

        if num_gen == 0:
            previous_energy_path = path_previous_energy_file['initial']
        else:
            previous_energy_path = path_previous_energy_file['not_initial']

        actual_energy_path = f'structure_files/generations/generation-{num_gen + 1:03d}/relaxed_structures/energy_ranking.txt'

        previous_energy = open(previous_energy_path, "r")
        previous_energy.readline()

        actual_energy = open(actual_energy_path, "r")
        actual_energy.readline()

        file_energy = open(f'structure_files/generations/generation-{num_gen + 1:03d}/final_structures/energy_ranking.txt', 'w')
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
                shutil.copy(f'structure_files/generations/generation-{num_gen + 1:03d}/relaxed_structures/' + struc_file, 
                                f'structure_files/generations/generation-{num_gen + 1:03d}/final_structures')
                    
                file_energy.write(line_actual_energy)
            else:
                shutil.copy(f'structure_files/generations/generation-{num_gen + 1:03d}/initial_structures/' + struc_file, 
                                f'structure_files/generations/generation-{num_gen + 1:03d}/final_structures')
                    
                file_energy.write(line_previous_energy)

        file_energy.close()        
        previous_energy.close()
        actual_energy.close()

        if print_terminal_outputs == True:
            gen_end_actual(num_gen + 1)

        num_gen = num_gen + 1

    if print_terminal_outputs == True:
        gen_end()

    os.mkdir('structure_files/generations/final_structures/')

    for struc_file in name_selected_structures:
        shutil.copy(f'structure_files/generations/generation-{num_generations:03d}/initial_structures/' + struc_file, 
                    'structure_files/generations/final_structures/' + struc_file)
            
    phases_and_energies = []
    final_energy_file = open(f'structure_files/generations/generation-{num_generations:03d}/final_structures/energy_ranking.txt', 'r')
    final_energy_file.readline()
    for it in range(len(name_selected_structures)):
        line_file = final_energy_file.readline()
        phases_and_energies_element = [line_file.split()[1], float(line_file.split()[2]), line_file.split()[3] + ' ' + line_file.split()[4]]
        phases_and_energies.append(phases_and_energies_element)
    final_energy_file.close()
        
    sorted_phases_and_energies = sorted(phases_and_energies, key=lambda x: x[1])

    final_energy_file = open('structure_files/generations/final_structures/energy_ranking.txt', 'w')
    if structure_file == 'poscar':
        final_energy_file.write('#       POSCAR-num       energy per atom (eV)       Space Group (Hermann-Mauguin)\n')
    elif structure_file == 'cif':
        final_energy_file.write('#       structure-num.cif       energy per atom (eV)       Space Group (Hermann-Mauguin)\n')

    num_phase = 1
    for phase in sorted_phases_and_energies:
        final_energy_file.write(f'{num_phase}       {phase[0]}       {phase[1]}       {phase[2]}\n')
        num_phase = num_phase + 1
        
    final_energy_file.close()

    deduplicate_file('structure_files/generations/final_structures/energy_ranking.txt', 'structure_files/generations/final_structures/', 
                     num_phase - 1, structure_file, prec_group_det)

    end_time = time.time()
    elapsed_time = end_time - start_time

    if print_terminal_outputs == True:
        print(f'Elapsed time: {elapsed_time:.3f} seconds')


    ### Close the log file
        
    if save_log_file == True:
        log_file_handle.close()
        sys.stdout = sys.__stdout__

    return


def change_exp_tuned_parameters(prominance, width, inputs_path):
    """
    Change the inputs_diffraction file with the new prominance and widht parameters

    Inputs:
        prominance: how much peak stands out respect the baseline
        width: width of the peaks
        inputs_path: path to inputs_diffraction file
    """

    string_array = []
    inputs = open(inputs_path, 'r')
    for num_lines in range(33):
        actual_line = inputs.readline()
        string_array.append(actual_line)
    inputs.close()

    inputs = open(inputs_path, 'w')
    for num_lines in range(26):
        inputs.write(string_array[num_lines])
    inputs.write(f'prominance_exp          =  {prominance}\n')
    inputs.write(f'width_exp               =  {width}\n')
    for num_lines in range(5):
        inputs.write(string_array[num_lines + 28])
    inputs.close()

    return