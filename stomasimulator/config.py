""" Handles application configuration """

from __future__ import print_function

import os

from .utils import utils as ut


class AppConfig(object):
    """ Store the application configuration """

    instance = None

    def __init__(self):
        if AppConfig.instance is None:
            AppConfig.instance = AppConfig._AppConfig()

    def __getattr__(self, item):
        return getattr(self.instance, item)

    class _AppConfig(object):
        def __init__(self):
            cfg_fn = self.get_config_filename()

            print('*' * 120)
            print("Reading configuration file '{}'...".format(cfg_fn))

            if not os.path.exists(cfg_fn):
                raise ValueError('The configuration file does not exist')

            if not os.path.isfile(cfg_fn):
                raise ValueError('The specified configuration file is not a file')

            config_dict = ut.read_yaml(cfg_fn)

            self.febio_exe = os.path.expanduser( config_dict['febio-exe'] )
            self.model_data_dir = os.path.join(self.get_root_dir(), config_dict['model-data-dir'])

            print('Finished reading the configuration file')
            print('')

            self.print_config()

            print('*' * 120)

        @staticmethod
        def get_root_dir():
            """ Full path to this directory """
            return os.path.dirname(os.path.abspath(__file__))

        def get_config_filename(self):
            """ The name of the application configuration file """
            return os.path.join(self.get_root_dir(), 'stomasimulator.yml')

        def print_config(self):
            """ Output the configuration """

            print('Configuration settings for stomasimulator...')
            print('FEBio executable     : {}'.format(self.febio_exe))
            print('Model data directory : {}'.format(self.model_data_dir))


APP_CONFIG = AppConfig()


class YAMLFileContents(object):
    """ Representation of the YAML file that stores model configuration """

    def __init__(self, file_name, file_contents):
        for key in ('observation', 'model', 'pressure', 'optimisation'):
            if key not in file_contents:
                raise KeyError('The YAML file does not contain the mandatory key, {}, at the top level'.format(key))

        self.file_name = file_name

        self.observations = file_contents['observation']
        self.models = file_contents['model']
        self.pressure = file_contents['pressure']
        self.optimisation = file_contents['optimisation']

    def __repr__(self):
        return '[filename={}, observations={}, models={}, pressure={}, optimisation={}'.format(
            self.file_name,
            self.observations.keys(),
            self.models.keys(),
            self.pressure,
            self.optimisation['keys'])


class YAMLFilesReader(object):
    """ Iterate over the YAML files in a directory """

    def __init__(self, dir_name):
        self.dir_name = dir_name
        self.file_names = ut.get_file_list(dir_name=self.dir_name,
                                           suffix_filters=('yml', 'yaml'))

    def __iter__(self):
        for file_name in self.file_names:
            full_file_name = os.path.join(self.dir_name, file_name)
            yaml_file_contents = ut.read_yaml(full_file_name)

            yield YAMLFileContents(file_name=file_name, file_contents=yaml_file_contents)


def get_file_contents_for_model_label(model_label):
    """
    Get the contents of the YAML file that contains the given model

    :param model_label: search for a model with this label
    :return: YAMLFileContents object
    :rtype: YAMLFileContents
    """

    yaml_files_reader = YAMLFilesReader(APP_CONFIG.model_data_dir)

    for yaml_file in yaml_files_reader:
        if model_label in yaml_file.models:
            return yaml_file

    raise KeyError('The model label, {}, was not found'.format(model_label))


def get_model_labels():
    """
    Get the model labels from the config file

    :return: a sorted list of model labels
    :rtype: list
    """
    model_labels = []

    yaml_files_reader = YAMLFilesReader(APP_CONFIG.model_data_dir)

    for yaml_file in yaml_files_reader:
        model_labels += yaml_file.models.keys()

    return sorted(model_labels)


def print_model_labels(filter_str):
    """
    Output the model labels from the config file
    :return:
    """

    yaml_files_reader = YAMLFilesReader(APP_CONFIG.model_data_dir)

    if filter_str:
        print("Filtering model labels using '{}'".format(filter_str))

    num_models = 0

    for yaml_file in yaml_files_reader:
        if filter_str:
            model_labels = [_ for _ in yaml_file.models.keys() if filter_str in _]
        else:
            model_labels = yaml_file.models.keys()

        model_labels.sort()

        if len(model_labels) > 0:
            print('Model labels for: {}'.format(yaml_file.file_name))
            for model_label in model_labels:
                print('   {}'.format(model_label))
                num_models += 1

    print('')
    print('Found {} model(s)'.format(num_models))
    print('')


if __name__ == '__main__':
    pass
