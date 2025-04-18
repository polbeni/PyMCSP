# Pol Benítez Colominas, Universitat Politècnica de Catalunya
# September 2023 - March 2025

# Diffraction module

import os
import csv
from copy import deepcopy
import multiprocessing

import numpy as np

from pymatgen.io.vasp import Poscar
from pymatgen.analysis.diffraction.xrd import XRDCalculator
from pymatgen.analysis.diffraction.neutron import NDCalculator

from functions_pymcsp import *
from terminal_outputs import *
from read_inputs import *


log_file = 'output_diffraction.log'

#### Load default inputs and look for new inputs in a YAML file

if os.path.exists('inputs.yaml'):
    config = InputsPyMCSP('inputs.yaml')
else:
    config = InputsPyMCSP()

path_structures = config.path_structures
structure_file = config.structure_file
name_exp_diff  = config.name_exp_diff
clean_substrate = config.clean_substrate
path_substrate = config.path_substrate
tolerance_subs = config.tolerance_subs
type_diffraction = config.type_diffraction
wl_xray = config.wl_xray
wl_neutron = config.wl_neutron
min_2theta = config.min_2theta
max_2theta = config.max_2theta
total_num_points = config.total_num_points
peak_width = config.peak_width
type_interpolation = config.type_interpolation
prominance_exp = config.prominance_exp
width_exp = config.width_exp
prop_vol = config.prop_vol
num_vols = config.num_vols
coef_quad = config.coef_quad
print_terminal_outputs = config.print_terminal_outputs
save_log_file = config.save_log_file


#### Log file

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

start_time = time.time()


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

    scale_vol = minimize_loss_vol_scale(structure, calculator, min_2theta, max_2theta, total_num_points, type_interpolation, peak_width, prop_vol, num_vols, exp_normalized_intensity)

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
        if print_terminal_outputs == True:
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

end_time = time.time()
elapsed_time = end_time - start_time

print(f'Elapsed time: {elapsed_time:.3f} seconds')

if save_log_file == True:
    log_file_handle.close()
    sys.stdout = sys.__stdout__