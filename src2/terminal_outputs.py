# Pol Benítez Colominas, Universitat Politècnica de Catalunya
# September 2023 - March 2024
# Version 1.0

# Functions for output messages in terminal

def initial_message():
    """
    Message to show when start the program execution
    """

    print('')
    print('o---------------------------------------------------o')
    print('|                      PyMCSP                       |')
    print('|           Crystal Structure Prediction            |')
    print('|               and Diffraction Study               |')
    print('|                                                   |')
    print('|                 Pol Benítez Colominas, 2023-2024  |')
    print('|             Universitat Politècnica de Catalunya  |')
    print('o---------------------------------------------------o')
    print('')

    return


def menu_display():
    """
    Displays the possible options of the program
    """

    print('')
    print('o---------------------------------------------------o')
    print('What do you want to do?')
    print('-Crystal Structure prediction:')
    print('   1: Perform crystal structure prediction')
    print('   2: Pressure computations for a given structures')
    print('   3: Generations loop for a given structures')
    print('-Diffraction study:')
    print('   4: Find the structural phase of experimental results')
    print('   5: Tune experimental diffraction curve conversion')
    print('-Others:')
    print('   9: About')
    print('   0: Exit')
    print('o---------------------------------------------------o')
    print('')


def option_1(inputs):
    """
    Text for option 1 (Perform crystal structure prediction)

    Inputs:
        inputs: array with all the necessary inputs from inputs_csp file
    """

    print('')
    print('o---------------------------------------------------o')
    print('You have selected: Perform crystal structure prediction')
    print('')
    print('The inputs are:')
    print(f'dimension               =  {int(inputs[0])}')
    print(f'atoms                   =  {inputs[1]}')
    print(f'stoichiometry           =  {inputs[2]}')
    print(f'num_max_atoms           =  {int(inputs[3])}')
    print(f'num_same_phase          =  {int(inputs[4])}')
    print(f'max_ionic_steps         =  {int(inputs[5])}')
    print(f'comp_pressure           =  {inputs[6]}')
    print(f'pressure                =  {float(inputs[7])}')
    print(f'num_volumes             =  {int(inputs[8])}')
    print(f'minimum_volume          =  {float(inputs[9])}')
    print(f'maximum_volume          =  {float(inputs[10])}')
    print(f'print_terminal_outputs  =  {inputs[11]}')
    print(f'save_log_file           =  {inputs[12]}')
    print(f'save_all_generations    =  {inputs[13]}')
    print(f'structure_file          =  {inputs[14]}')
    print(f'prec_group_det          =  {float(inputs[15])}')
    print(f'num_generations         =  {int(inputs[16])}')
    print(f'surviving_phases        =  {float(inputs[17])}')
    print(f'max_disp                =  {float(inputs[18])}')
    print('')
    print('Are they correct? (write yes or no)')
    print('o---------------------------------------------------o')
    print('')


def option_2(inputs):
    """
    Text for option 2 (Pressure computations for a given structures)

    Inputs:
        inputs: array with all the necessary inputs from inputs_csp file
    """

    print('')
    print('o---------------------------------------------------o')
    print('You have selected: Pressure computations for a given structures')
    print('')
    print('The inputs are:')
    print(f'max_ionic_steps         =  {int(inputs[5])}')
    print(f'comp_pressure           =  {inputs[6]}')
    print(f'pressure                =  {float(inputs[7])}')
    print(f'num_volumes             =  {int(inputs[8])}')
    print(f'minimum_volume          =  {float(inputs[9])}')
    print(f'maximum_volume          =  {float(inputs[10])}')
    print(f'print_terminal_outputs  =  {inputs[11]}')
    print(f'save_log_file           =  {inputs[12]}')
    print(f'save_all_generations    =  {inputs[13]}')
    print(f'structure_file          =  {inputs[14]}')
    print(f'prec_group_det          =  {float(inputs[15])}')
    print('')
    print('Are they correct? (write yes or no)')
    print('o---------------------------------------------------o')
    print('')

def option_3(inputs):
    """
    Text for option 3 (Generations loop for a given structures)

    Inputs:
        inputs: array with all the necessary inputs from inputs_csp file
    """

    print('')
    print('o---------------------------------------------------o')
    print('You have selected: Generations loop for a given structures')
    print('')
    print('The inputs are:')
    print(f'max_ionic_steps         =  {int(inputs[5])}')
    print(f'print_terminal_outputs  =  {inputs[11]}')
    print(f'save_log_file           =  {inputs[12]}')
    print(f'save_all_generations    =  {inputs[13]}')
    print(f'structure_file          =  {inputs[14]}')
    print(f'prec_group_det          =  {float(inputs[15])}')
    print(f'num_generations         =  {int(inputs[16])}')
    print(f'surviving_phases        =  {float(inputs[17])}')
    print(f'max_disp                =  {float(inputs[18])}')
    print('')
    print('Are they correct? (write yes or no)')
    print('o---------------------------------------------------o')
    print('')


def option_4(inputs):
    """
    Text for option 4 (Find the structural phase of experimental results)

    Inputs:
        inputs: array with all the necessary inputs from inputs_csp file
    """

    print('')
    print('o---------------------------------------------------o')
    print('You have selected: Find the structural phase of experimental results')
    print('')
    print('The inputs are:')
    print(f'path_structures         =  {inputs[0]}')
    print(f'structure_file          =  {inputs[1]}')
    print(f'name_exp_diff           =  {inputs[2]}')
    print(f'clean_substrate         =  {inputs[3]}')
    print(f'path_substrate          =  {inputs[4]}')
    print(f'tolerance_subs          =  {float(inputs[5])}')
    print(f'type_diffraction        =  {inputs[6]}')
    print(f'wl_xray                 =  {inputs[7]}')
    print(f'wl_neutron              =  {float(inputs[8])}')
    print(f'min_2theta              =  {float(inputs[9])}')
    print(f'max_2theta              =  {float(inputs[10])}')
    print(f'total_num_points        =  {int(inputs[11])}')
    print(f'peak_width              =  {float(inputs[12])}')
    print(f'type_interpolation      =  {inputs[13]}')
    print(f'prominance_exp          =  {float(inputs[14])}')
    print(f'width_exp               =  {float(inputs[15])}')
    print(f'prop_vol                =  {float(inputs[16])}')
    print(f'num_vols                =  {int(inputs[17])}')
    print(f'coef_quad               =  {float(inputs[18])}')
    print('')
    print('Are they correct? (write yes or no)')
    print('o---------------------------------------------------o')
    print('')


def option_5():
    """
    Text for option 5 (Tune experimental diffraction curve conversion)
    """

    print('')
    print('o---------------------------------------------------o')
    print('You have selected: Tune experimental diffraction curve conversion')
    print('o---------------------------------------------------o')
    print('')


def option_5_parameters():
    """
    Get the parameters to tune the curve
    """

    print('')
    print('Select a value for prominance (default: 100):')
    prominance = input()
    print('Select a value for width (default: 5):')
    width = input()
    print('')

    return float(prominance), float(width)


def option_5_final():
    """
    Final text for option 5
    """

    print('')
    print('o---------------------------------------------------o')
    print('The values of prominance_exp and width_exp have been changed successfully!')
    print('o---------------------------------------------------o')
    print('')


def about():
    """
    Text for about option
    """

    print('')
    print('o---------------------------------------------------o')
    print('This code has been developed in Universitat Politècnica')
    print('de Catalunya by:')
    print('   Pol Benítez Colominas (pol.benitez@upc.edu)')
    print('   Claudio Cazorla Silva (claudio.cazorla@upc.edu)')
    print('If you have any comment or bug you would like to report,')
    print('please do not hesitate to contact us!')
    print('')
    print('o---------------------------------------------------o')
    print('Press any key to finish')
    print('o---------------------------------------------------o')
    print('')


def error_message():
    """
    Message if the provided input is not valid
    """

    print('')
    print('o---------------------------------------------------o')
    print('The provided input is not valid, please try again')
    print('o---------------------------------------------------o')
    print('')


def final_message():
    """
    Message to show when finish the program execution
    """

    print('')
    print('o---------------------------------------------------o')
    print('|                   See you soon!                   |')
    print('o---------------------------------------------------o')
    print('')

    return


def simulation_information(atoms, stoichiometry, num_max, gen):
    """
    Provides information about the simulation

    Inptus:
        atoms: array with the atoms of the material
        stoichiometry: nmber of atoms for each element
        num_max: maximum number of atoms in the unit cell
        gen: number of generations
    """

    material = ''

    for x in range(len(atoms)):
        material = material + str(atoms[x]) + str(stoichiometry[x])

    print('')
    print('The simulation has started!')
    print(f'Looking for phases of the material {material}.')
    print(f'The maximum number of atoms per unit cell is {num_max}.')
    print(f'The number of generations is {gen}.')
    print('')

    return


def struc_gen_ini():
    """
    Starting the structure generation
    """

    print('')
    print('Starting the structure generation.')
    print('')

    return

    
def struc_gen_end():
    """
    End the structure generation
    """

    print('')
    print('Structures generated!')
    print('')

    return


def relax_ini():
    """
    Starting the structure relaxation
    """

    print('')
    print('Relaxing the structures.')
    print('')

    return

    
def relax_end():
    """
    End the structure relaxation
    """

    print('')
    print('Structures relaxed!')
    print('')

    return


def pressure_ini():
    """
    Starting the pressure calculations
    """

    print('')
    print('Determining volumes that minimize the enthalpy for pressure study.')
    print('')

    return

    
def pressure_end():
    """
    End the pressure calculations
    """

    print('')
    print('Energies determined for the structures with pressure!')
    print('')

    return


def gen_ini():
    """
    Starting generations
    """

    print('')
    print('Entering the generations loop.')
    print('')

    return

    
def gen_end():
    """
    End generations
    """

    print('')
    print('Generations loop finshed!')
    print('')

    return


def gen_ini_actual(gen):
    """
    Starting generations

    Inputs:
        gen: number of the generation
    """

    print('')
    print(f'Entering the generation number {gen}.')
    print('')

    return

    
def gen_end_actual(gen):
    """
    End generations

    Inputs:
        gen: number of the generation
    """

    print('')
    print(f'Generation number {gen} finshed!')
    print('')

    return


def finished_execution():
    """
    Message to show when the program execution finishes
    """

    print('')
    print('o---------------------------------------------------o')
    print('|            The simulation has finished!           |')
    print('o---------------------------------------------------o')
    print('')

    return