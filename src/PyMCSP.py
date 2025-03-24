# Pol Benítez Colominas, Universitat Politècnica de Catalunya
# September 2023 - March 2025

# Main script

import os
import sys
import shutil
import time
import subprocess

from functions_pymcsp import *
from terminal_outputs import *

import warnings

warnings.simplefilter("ignore") # To suppress warnings for clearer output

initial_message()

finish_condition = False

while finish_condition == False:
    menu_display()

    user_input = input()

    if user_input == '1':
        inputs_var = read_variables('inputs.yaml')
        option_1(inputs_var)
        correct = input()
        if correct == 'yes':
            csp_study(inputs_var)
            finish_condition = True
        elif correct == 'no':
            finish_condition = True
        else:
            error_message()

    elif user_input == '2':
        inputs_var = read_variables('inputs.yaml')
        option_2(inputs_var)
        correct = input()
        if correct == 'yes':
            pressure_computations(inputs_var)
            finish_condition = True
        elif correct == 'no':
            finish_condition = True
        else:
            error_message()

    elif user_input == '3':
        inputs_var = read_variables('inputs.yaml')
        option_3(inputs_var)
        correct = input()
        if correct == 'yes':
            generations_loop(inputs_var)
            finish_condition = True
        elif correct == 'no':
            finish_condition = True
        else:
            error_message()

    elif user_input == '4':
        inputs_var = read_variables('inputs.yaml')
        option_4(inputs_var)
        correct = input()
        if correct == 'yes':
            subprocess.run(['python3', 'diffraction_module.py'])
            finish_condition = True
        elif correct == 'no':
            finish_condition = True
        else:
            error_message()

    elif user_input == '5':
        inputs_var = read_variables('inputs.yaml')
        option_5()
        correct_parameters = False
        while correct_parameters == False:
            prominance, width = option_5_parameters()
            plot_exp_tuned(prominance, width, inputs_var[2])
            print('')
            print('Are the results okay? (answer yes or no)')
            agreement = input()
            if agreement == 'yes':
                correct_parameters = True
                finish_condition = True
                option_5_final()
            elif agreement == 'no':
                correct_parameters = False
            else:
                error_message()
                correct_parameters = False
            print('')
            

    elif user_input == '9':
        about()
        input()
        finish_condition = True

    elif user_input == '0':
        finish_condition = True

    else: 
        error_message()

final_message()