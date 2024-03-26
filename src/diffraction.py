# Pol Benítez Colominas, Universitat Politècnica de Catalunya
# September 2023 - March 2024
# Version 0.3

# Code to perform diffraction study


#### Import libraries

import csv
from copy import deepcopy

import numpy as np
from scipy.interpolate import interp1d
from scipy import signal
from scipy import integrate

from pymatgen.io.vasp import Poscar
from pymatgen.analysis.diffraction.xrd import XRDCalculator
from pymatgen.analysis.diffraction.neutron import NDCalculator

import matplotlib.pyplot as plt

##### IMPORTANT #####
# XRDCalculator and NDCalculator have wavelength variable, for now using the default values
# read https://pymatgen.org/pymatgen.analysis.diffraction.html

structure = Poscar.from_file('structure_files/initial_structures/relaxed_structures/POSCAR-2223').structure

##### Variables inputs #####
path_structures = 'structure_files/initial_structures/relaxed_structures/'
structure_file = 'poscar'
name_exp_diff = 'diff_exp.csv' # the experimental diffractogram has to be a csv file with two columns, 2theta and intensity
type_diffraction = 'x-ray'
wl_xray = 'CuKa' # default value
wl_neutron = 1.54184 # default value
min_2theta = 20
max_2theta = 60
total_num_points = 2000 # recommended not decrease the value
peak_width = 0.25 
type_interpolation = 'linear' # recommended not touch
prominance_exp = 100
width_exp = 5
displacement_range = 5
############################

def curve_from_peaks(peaks_x, peaks_y, min_x, max_x, num_points, kind_interpolation):
    """
    Returns a curve from a given set of peaks

    Inputs:
        peaks_x: value of the peaks in x range
        peaks_y: value of the peaks in y range
        min_x: minimum value of the x range of the interpolation
        max_x: maximum value of the x range of the interpolation
        num_points: number of points of the interpolation
        kind_interpolation: type of interpolation performed
    """

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


def compute_loss_factor(range_x, intensity_theoretical, intesity_experimental, num_peaks_theo, num_peaks_exp):
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
    """

    function_to_integrate = np.zeros(len(range_x))
    function_to_integrate[:] = np.abs(intesity_experimental[:] - intensity_theoretical[:])**2

    T_result = integrate.trapezoid(function_to_integrate, range_x) + np.abs(num_peaks_theo - num_peaks_exp)**2

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


def minimize_loss_vol_scale(structure, min_x, max_x, num_points, kind_interpolation, prop_vol, num_vols, intesity_experimental):
    """
    Scales the volume of the structure in the desired interval and finds the structure with minimum loss factor

    Inputs:
        structure: structure pymatgen object to study
        min_x: minimum value of the x range of the interpolation
        max_x: maximum value of the x range of the interpolation
        num_points: number of points of the interpolation
        kind_interpolation: type of interpolation performed
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

        twotheta, intensity = curve_from_peaks(two_theta_peaks, intensity_peaks, min_x, max_x, num_points, kind_interpolation)

        normalized_intensity = normalize_curve(twotheta, intensity)

        loss_factor = compute_loss_factor_minim(twotheta, normalized_intensity, intesity_experimental)
        loss_array.append(loss_factor)

    min_loss_index = loss_array.index(min(loss_array))
    vol_to_minimize = vol_array[min_loss_index]

    return vol_to_minimize

    


# theoretical

if type_diffraction == 'x-ray':
    calculator = XRDCalculator(wavelength=wl_xray)
elif type_diffraction == 'neutron':
    calculator = NDCalculator(wavelength=wl_neutron)
else:
    print('Not valid type_diffraction value')


pattern = calculator.get_pattern(structure)

two_theta_peaks = pattern.x
intensity_peaks = pattern.y

twotheta, intensity = curve_from_peaks(two_theta_peaks, intensity_peaks, min_2theta, max_2theta, total_num_points, type_interpolation)

normalized_intensity = normalize_curve(twotheta, intensity)


# experimental

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

exp_twotheta, exp_intensity = curve_from_peaks(exp_peaks_x, exp_peaks_y, min_2theta, max_2theta, total_num_points, type_interpolation)

exp_normalized_intensity = normalize_curve(exp_twotheta, exp_intensity)


#loss_factor = compute_loss_factor(exp_twotheta, normalized_intensity, exp_normalized_intensity)

#print(loss_factor)

scale_vol = minimize_loss_vol_scale(structure, min_2theta, max_2theta, total_num_points, type_interpolation, 0.05, 20, exp_normalized_intensity)

distorted = deepcopy(structure)
distorted.scale_lattice(distorted.volume * scale_vol)

pattern = calculator.get_pattern(distorted)
two_theta_peaks = pattern.x
intensity_peaks = pattern.y

twotheta, intensity = curve_from_peaks(two_theta_peaks, intensity_peaks, min_2theta, max_2theta, total_num_points, type_interpolation)

normalized_intensity = normalize_curve(twotheta, intensity)

num_peaks_theo = sum(1 for x in two_theta_peaks if min_2theta <= x <= max_2theta)
num_peaks_exp = sum(1 for x in exp_peaks_x if min_2theta <= x <= max_2theta)

loss_factor = compute_loss_factor(exp_twotheta, normalized_intensity, exp_normalized_intensity, num_peaks_theo, num_peaks_exp)
print(num_peaks_theo, num_peaks_exp)
print(loss_factor)

#print(loss_factor)

#plt.figure()
#plt.plot(twotheta, normalized_intensity)
#plt.plot(exp_twotheta, exp_normalized_intensity)
#plt.show()