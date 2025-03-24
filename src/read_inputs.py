# Pol Benítez Colominas, Universitat Politècnica de Catalunya
# March 2025 - March 2025

# Default inputs and read inputs

import yaml

class InputsPyMCSP:
    """
    PyMCSP inputs, it uses the default values or in the inputs.yaml file
    """

    # Define default values for the inputs
    DEFAULT_CONFIG = {
        'dimension': 3,
        'atoms': ['Ga', 'As'],
        'stoichiometry': [1, 1],
        'num_max_atoms': 20,
        'num_same_phase': 10,
        'max_ionic_steps': 250,
        'pressure': 1e8,
        'num_volumes': 30,
        'minimum_volume': 0.8,
        'maximum_volume': 1.2,
        'print_terminal_outputs': True,
        'save_log_file': True,
        'structure_file': 'poscar',
        'prec_group_det': 1e-5,
        'num_generations': 3,
        'surviving_phases': 0.2,
        'max_disp': 0.1,
        'struc_path': 'structure_files/relaxed_structures',
        'retrain': False,
        'retrain_path': 'path-retreined',
        'restricted_phases': False,
        'restricted_list': [1, 22, 47, 124],
        'path_structures': 'structure_files/relaxed_structures/',
        'name_exp_diff': 'diff_exp.csv',
        'clean_substrate': False,
        'path_substrate': 'SnO2.poscar',
        'tolerance_subs': 0.1,
        'type_diffraction': 'x-ray',
        'wl_xray': 'CuKa',
        'wl_neutron': 1.54184,
        'min_2theta': 20,
        'max_2theta': 60,
        'total_num_points': 2000,
        'peak_width': 0.25,
        'type_interpolation': 'linear',
        'prominance_exp': '100.0',
        'width_exp': 5.0,
        'prop_vol': 0.05,
        'num_vols': 20,
        'coef_quad': 0.01
    }

    def __init__(self, config_file: str = None):
        """
        Initialize with default values, optionally loading from YAML file
        """

        # Copy the default values
        self.config = self.DEFAULT_CONFIG.copy()

        # Update with file contents if provided
        if config_file:
            self.load_config(config_file)

    def load_config(self, file_path: str):
        """
        Load YAML file if provided and update specified values
        """

        with open(file_path, 'r') as f:
            user_config = yaml.safe_load(f) or {}
            self.config.update(user_config)

    def __getattr__(self, name):
        """
        Allow attribute-style access
        """
        
        try:
            return self.config[name]
        except KeyError as e:
            raise AttributeError(f"No such config option: {name}") from e