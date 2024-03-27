# Pol Benítez Colominas, Universitat Politècnica de Catalunya
# September 2023 - March 2024
# Version 0.4

# Main script to perfrom diffraction study


#### Import libraries

import os
import csv
from copy import deepcopy
import multiprocessing
from functools import partial

import numpy as np
from scipy.interpolate import interp1d
from scipy import signal
from scipy import integrate

from pymatgen.io.vasp import Poscar
from pymatgen.analysis.diffraction.xrd import XRDCalculator
from pymatgen.analysis.diffraction.neutron import NDCalculator


#### Read inputs from inputs file

path_structures = None
structure_file = None
name_exp_diff  = None
clean_substrate = None
path_substrate = None
tolerance_subs = None
type_diffraction = None
wl_xray = None
wl_neutron = None
min_2theta = None
max_2theta = None
total_num_points = None
peak_width = None
type_interpolation = None
prominance_exp = None
width_exp = None
prop_vol = None
num_vols = None
coef_quad = None

inputs = open('inputs_diffraction', "r")

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

path_structures = variables[0]
structure_file = variables[1]
name_exp_diff  = variables[2]
clean_substrate = variables[3]
path_substrate = variables[4]
tolerance_subs = float(variables[5])
type_diffraction = variables[6]
wl_xray = variables[7]
wl_neutron = float(variables[8])
min_2theta = float(variables[9])
max_2theta = float(variables[10])
total_num_points = int(variables[11])
peak_width = float(variables[12])
type_interpolation = variables[13]
prominance_exp = float(variables[14])
width_exp = float(variables[15])
prop_vol = float(variables[16])
num_vols = int(variables[17])
coef_quad = float(variables[18])


#### Functions

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
    substrate = Poscar.from_file(path_substrate).structure

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

structures_list = []
final_energy_file = open(f'{path_structures}/energy_ranking.txt', 'r')
final_energy_file.readline()
for it in range(num_structures):
    line_file = final_energy_file.readline()
    structure_element = [line_file.split()[1], None, line_file.split()[3] + ' ' + line_file.split()[4]]
    structures_list.append(structure_element)
final_energy_file.close()

diffraction_results = open('diffraction_results.txt', 'w')
if structure_file == 'poscar':
    diffraction_results.write('#       POSCAR-num       Loss factor       Space Group (Hermann-Mauguin)\n')
elif structure_file == 'cif':
    diffraction_results.write('#       structure-num.cif       Loss factor       Space Group (Hermann-Mauguin)\n')


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
        print(f'Loss factor computed for the structure number {counter.value} of {num_structures} structures.')

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
    diffraction_results.write(f'{num_phase}       {phase[0]}       {phase[1]}       {phase[2]}\n')
    num_phase = num_phase + 1

diffraction_results.close()