# Pol Benítez Colominas, Universitat Politècnica de Catalunya
# September 2023 - September 2023
# Version 0.2

# Functions for output messages in terminal file

def initial_message():
    """
    Message to show when start the program execution
    """

    print('')
    print('o---------------------------------------------------o')
    print('|                      PyMCSP                       |')
    print('|                                                   |')
    print('|                      Pol Benítez Colominas, 2023  |')
    print('|             Universitat Politècnica de Catalunya  |')
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


def final_message():
    """
    Message to show when the program execution finishes
    """

    print('')
    print('o---------------------------------------------------o')
    print('|            The simulation has finished!           |')
    print('o---------------------------------------------------o')
    print('')

    return