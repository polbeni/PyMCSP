# Pol Benítez Colominas, Universitat Politècnica de Catalunya
# September 2023 - March 2025

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
    print('|                 Pol Benítez Colominas, 2023-2025  |')
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
    print(f'dimension               =  {inputs.dimension}')
    print(f'atoms                   =  {inputs.atoms}')
    print(f'stoichiometry           =  {inputs.stoichiometry}')
    print(f'num_max_atoms           =  {inputs.num_max_atoms}')
    print(f'num_same_phase          =  {inputs.num_same_phase}')
    print(f'max_ionic_steps         =  {inputs.max_ionic_steps}')
    print(f'print_terminal_outputs  =  {inputs.print_terminal_outputs}')
    print(f'save_log_file           =  {inputs.save_log_file}')
    print(f'structure_file          =  {inputs.structure_file}')
    print(f'prec_group_det          =  {inputs.prec_group_det}')
    print(f'retrain                 =  {inputs.retrain}')
    if inputs.retrain == True:
        print(f'retrain_path            =  {inputs.retrain_path}')
    print(f'restricted_phases       =  {inputs.restricted_phases}')
    if inputs.restricted_phases == True:
        print(f'restricted_list         =  {inputs.restricted_list}')
    print(f'model                   =  {inputs.model}')
    if inputs.model == 'MACE':
        print(f'mace_model              =  {inputs.mace_model}')
        print(f'device_mace             =  {inputs.device_mace}')
        print(f'fmax_mace               =  {inputs.fmax_mace}')
        print(f'precision_relax_mace    =  {inputs.precision_relax_mace}')
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
    print(f'pressure                =  {inputs.pressure}')
    print(f'num_volumes             =  {inputs.num_volumes}')
    print(f'minimum_volume          =  {inputs.minimum_volume}')
    print(f'maximum_volume          =  {inputs.maximum_volume}')
    print(f'print_terminal_outputs  =  {inputs.print_terminal_outputs}')
    print(f'save_log_file           =  {inputs.save_log_file}')
    print(f'structure_file          =  {inputs.structure_file}')
    print(f'prec_group_det          =  {inputs.prec_group_det}')
    print(f'struc_path              =  {inputs.struc_path}')
    print(f'retrain                 =  {inputs.retrain}')
    if inputs.retrain == True:
        print(f'retrain_path            =  {inputs.retrain_path}')
    print(f'model                   =  {inputs.model}')
    if inputs.model == 'MACE':
        print(f'mace_model              =  {inputs.mace_model}')
        print(f'device_mace             =  {inputs.device_mace}')
        print(f'fmax_mace               =  {inputs.fmax_mace}')
        print(f'precision_relax_mace    =  {inputs.precision_relax_mace}')
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
    print(f'max_ionic_steps         =  {inputs.max_ionic_steps}')
    print(f'print_terminal_outputs  =  {inputs.print_terminal_outputs}')
    print(f'save_log_file           =  {inputs.save_log_file}')
    print(f'structure_file          =  {inputs.structure_file}')
    print(f'prec_group_det          =  {inputs.prec_group_det}')
    print(f'num_generations         =  {inputs.num_generations}')
    print(f'surviving_phases        =  {inputs.surviving_phases}')
    print(f'max_disp                =  {inputs.max_disp}')
    print(f'struc_path              =  {inputs.struc_path}')
    print(f'retrain                 =  {inputs.retrain}')
    if inputs.retrain == True:
        print(f'retrain_path            =  {inputs.retrain_path}')
    print(f'model                   =  {inputs.model}')
    if inputs.model == 'MACE':
        print(f'mace_model              =  {inputs.mace_model}')
        print(f'device_mace             =  {inputs.device_mace}')
        print(f'fmax_mace               =  {inputs.fmax_mace}')
        print(f'precision_relax_mace    =  {inputs.precision_relax_mace}')
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
    print(f'path_structures         =  {inputs.path_structures}')
    print(f'structure_file          =  {inputs.structure_file}')
    print(f'name_exp_diff           =  {inputs.name_exp_diff}')
    print(f'clean_substrate         =  {inputs.clean_substrate}')
    if inputs.clean_substrate == True:
        print(f'path_substrate          =  {inputs.path_substrate}')
        print(f'tolerance_subs          =  {inputs.tolerance_subs}')
    print(f'type_diffraction        =  {inputs.type_diffraction}')
    print(f'wl_xray                 =  {inputs.wl_xray}')
    print(f'wl_neutron              =  {inputs.wl_neutron}')
    print(f'min_2theta              =  {inputs.min_2theta}')
    print(f'max_2theta              =  {inputs.max_2theta}')
    print(f'total_num_points        =  {inputs.total_num_points}')
    print(f'peak_width              =  {inputs.peak_width}')
    print(f'type_interpolation      =  {inputs.type_interpolation}')
    print(f'prominance_exp          =  {inputs.prominance_exp}')
    print(f'width_exp               =  {inputs.width_exp}')
    print(f'prop_vol                =  {inputs.prop_vol}')
    print(f'num_vols                =  {inputs.num_vols}')
    print(f'coef_quad               =  {inputs.coef_quad}')
    print(f'print_terminal_outputs  =  {inputs.print_terminal_outputs}')
    print(f'save_log_file           =  {inputs.save_log_file}')
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


def simulation_information(atoms, stoichiometry, num_max):
    """
    Provides information about the simulation

    Inptus:
        atoms: array with the atoms of the material
        stoichiometry: nmber of atoms for each element
        num_max: maximum number of atoms in the unit cell
    """

    material = ''

    for x in range(len(atoms)):
        material = material + str(atoms[x]) + str(stoichiometry[x])

    print('')
    print('The simulation has started!')
    print(f'Looking for phases of the material {material}.')
    print(f'The maximum number of atoms per unit cell is {num_max}.')
    print('')

    return


def simulation_information_pressure(pressure):
    """
    Provides information about the pressure computations

    Inptus:
        pressure: pressure in Pa for pressure computations
    """

    print('')
    print('The simulation has started!')
    print(f'The selected pressure is {pressure} Pa.')
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


def relax_struc(struc):
    """
    Relaxing the structure number #
    """

    print('')
    print(f'Relaxing the structure: {struc}')
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


def press_ini():
    """
    Starting the pressure computation
    """

    print('')
    print('Starting the pressure computation.')
    print('')

    return


def press_struc(struc):
    """
    Computing pressure of the structure number #
    """

    print('')
    print(f'Computing pressure of the structure: {struc}')
    print('')

    return


def press_end():
    """
    End the pressure computation
    """

    print('')
    print('Structures with pressure determined!')
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