# Pol Benítez Colominas, Universitat Politècnica de Catalunya
# September 2023 - March 2024
# Version 1.0

# Functions file

import os
import sys
import shutil
import time
from copy import deepcopy
import csv
import multiprocessing
from functools import partial

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy import signal
from scipy import integrate

from pymatgen.core import Structure
from pymatgen.io.vasp import Poscar
from pymatgen.io.cif import CifParser
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.diffraction.xrd import XRDCalculator
from pymatgen.analysis.diffraction.neutron import NDCalculator

from pyxtal import pyxtal

from m3gnet.models import Relaxer

import warnings

for category in (UserWarning, DeprecationWarning):
    warnings.filterwarnings("ignore", category=category, module="tensorflow")


def read_variables_csp(path_file):
    """
    Read the variables from the inputs_csp file and returns an array with their values

    Inputs:
        path_file: path of the inputs file
    """

    inputs = open(path_file, "r")

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

    return variables


def read_variables_diffraction(path_file):
    """
    Read the variables from the inputs_diffraction file and returns an array with their values

    Inputs:
        path_file: path of the inputs file
    """

    inputs = open(path_file, "r")

    variables = [None]*19

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


def minimize_loss_vol_scale(structure, min_x, max_x, num_points, kind_interpolation, peak_width, prop_vol, num_vols, intesity_experimental):
    """
    Scales the volume of the structure in the desired interval and finds the structure with minimum loss factor

    Inputs:
        structure: structure pymatgen object to study
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
    for num_lines in range(31):
        actual_line = inputs.readline()
        string_array.append(actual_line)
    inputs.close()

    inputs = open(inputs_path, 'w')
    for num_lines in range(26):
        inputs.write(string_array[num_lines])
    inputs.write(f'prominance_exp          =  {prominance}\n')
    inputs.write(f'width_exp               =  {width}\n')
    for num_lines in range(3):
        inputs.write(string_array[num_lines + 28])
    inputs.close()

    return


def diffraction_study(inputs):
    """
    Performs difraction study

    Inputs:
        inputs: array with all the input variables from inputs_diffraction file
    """

    path_structures = inputs[0]
    structure_file = inputs[1]
    name_exp_diff  = inputs[2]
    clean_substrate = inputs[3]
    path_substrate = inputs[4]
    tolerance_subs = float(inputs[5])
    type_diffraction = inputs[6]
    wl_xray = inputs[7]
    wl_neutron = float(inputs[8])
    min_2theta = float(inputs[9])
    max_2theta = float(inputs[10])
    total_num_points = int(inputs[11])
    peak_width = float(inputs[12])
    type_interpolation = inputs[13]
    prominance_exp = float(inputs[14])
    width_exp = float(inputs[15])
    prop_vol = float(inputs[16])
    num_vols = int(inputs[17])
    coef_quad = float(inputs[18])

    #### Load calculator

    if type_diffraction == 'x-ray':
        calculator = XRDCalculator(wavelength=wl_xray)
    elif type_diffraction == 'neutron':
        calculator = NDCalculator(wavelength=wl_neutron)
    else:
        print('Not valid type_diffraction value')


    #### Experimental diffractogram

    exp_2theta = []
    exp_int = []

    with open(name_exp_diff, 'r') as file:
        reader = csv.reader(file)

        has_header = csv.Sniffer().has_header(file.read(1024))
        file.seek(0)
        
        if has_header:
            next(reader)  # Skip the header

        for row in reader:
            exp_2theta.append(float(row[0]))  
            exp_int.append(float(row[1]))  

    exp_2theta = np.array(exp_2theta)
    exp_int = np.array(exp_int)

    exp_peaks_x, exp_peaks_y = peaks_experimental(exp_2theta, exp_int, prominance_exp, width_exp)

    if clean_substrate == True:
        if structure_file == 'poscar':
            substrate = Poscar.from_file(path_substrate).structure
        elif structure_file == 'cif':
            parser = CifParser(path_substrate)
            substrate = parser.get_structures()[0]

        pattern = calculator.get_pattern(substrate)

        two_theta_peaks = pattern.x

        clean_x, clean_y = remove_substrate(exp_peaks_x, exp_peaks_y, two_theta_peaks, tolerance_subs)

        exp_twotheta, exp_intensity = curve_from_peaks(clean_x, clean_y, min_2theta, max_2theta, total_num_points, type_interpolation, peak_width)
    else:
        exp_twotheta, exp_intensity = curve_from_peaks(exp_peaks_x, exp_peaks_y, min_2theta, max_2theta, total_num_points, type_interpolation, peak_width)

    exp_normalized_intensity = normalize_curve(exp_twotheta, exp_intensity)
    

    #### Theoretical loop
        
    num_structures = 0

    for path in os.listdir(path_structures):
        if os.path.isfile(os.path.join(path_structures, path)):
            num_structures = num_structures + 1

    num_structures = num_structures - 1
    print(num_structures)

    structures_list = []
    final_energy_file = open(f'{path_structures}/energy_ranking.txt', 'r')
    final_energy_file.readline()
    for it in range(num_structures):
        line_file = final_energy_file.readline()
        structure_element = [line_file.split()[1], None, float(line_file.split()[2]), line_file.split()[3] + ' ' + line_file.split()[4]]
        structures_list.append(structure_element)
    final_energy_file.close()

    diffraction_results = open('diffraction_results.txt', 'w')
    if structure_file == 'poscar':
        diffraction_results.write('#       POSCAR-num       Loss factor      energy per atom (eV)       Space Group (Hermann-Mauguin)\n')
    elif structure_file == 'cif':
        diffraction_results.write('#       structure-num.cif       Loss factor      energy per atom (eV)       Space Group (Hermann-Mauguin)\n')


    def process_item(struc):
        structure = Poscar.from_file(f'{path_structures}/{struc[0]}').structure

        pattern = calculator.get_pattern(structure)

        two_theta_peaks = pattern.x
        intensity_peaks = pattern.y

        twotheta, intensity = curve_from_peaks(two_theta_peaks, intensity_peaks, min_2theta, max_2theta, total_num_points, type_interpolation, peak_width)

        normalized_intensity = normalize_curve(twotheta, intensity)

        scale_vol = minimize_loss_vol_scale(structure, min_2theta, max_2theta, total_num_points, type_interpolation, peak_width, prop_vol, num_vols, exp_normalized_intensity)

        distorted = deepcopy(structure)
        distorted.scale_lattice(distorted.volume * scale_vol)

        pattern = calculator.get_pattern(distorted)
        two_theta_peaks = pattern.x
        intensity_peaks = pattern.y

        peaks_x = []
        peaks_y = []
        for peak in range(len(intensity_peaks)):
            if intensity_peaks[peak] >= 1e-1:
                peaks_x.append(two_theta_peaks[peak])
                peaks_y.append(intensity_peaks[peak])

        twotheta, intensity = curve_from_peaks(two_theta_peaks, intensity_peaks, min_2theta, max_2theta, total_num_points, type_interpolation, peak_width)

        normalized_intensity = normalize_curve(twotheta, intensity)

        num_peaks_theo = sum(1 for x in peaks_x if min_2theta <= x <= max_2theta)
        num_peaks_exp = sum(1 for x in peaks_y if min_2theta <= x <= max_2theta)

        loss_factor = compute_loss_factor(exp_twotheta, normalized_intensity, exp_normalized_intensity, num_peaks_theo, num_peaks_exp, coef_quad)

        struc[1] = loss_factor

        with counter.get_lock():
            counter.value += 1
            print(f'Loss factor computed for the structure number {counter.value} of {num_structures}.')

        return struc

    def init_counter(counter_):
        global counter
        counter = counter_

    def get_counter_lock():
        return counter.get_lock()

    if __name__ == "__main__":
        elements = structures_list

    num_iteration = 1
    with multiprocessing.Pool(initializer=init_counter, initargs=(multiprocessing.Value('i', 0),)) as pool:
        results = pool.map(process_item, elements)

        for i, result in enumerate(results):
            elements[i] = result
        
    sorted_structures = sorted(elements, key=lambda x: x[1])

    num_phase = 1
    for phase in sorted_structures:
        diffraction_results.write(f'{num_phase}       {phase[0]}       {phase[1]}       {phase[2]}       {phase[3]}\n')
        num_phase = num_phase + 1

    diffraction_results.close()

    return