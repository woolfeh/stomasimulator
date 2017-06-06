"""
Utility functions
"""

from __future__ import print_function

import datetime
import os
import sys

# pylint: disable=import-error
import yaml


def module_exists(module_name):
    """
    Does the specified module exist?
    :param module_name:
    :return: boolean
    """
    try:
        __import__(module_name)
    except ImportError:
        return False
    else:
        return True


def get_datetime_now():
    """
    get time now
    :return:
    """
    return datetime.datetime.now()


def get_elapsed_millis(start_dt):
    """
    Calculate elapsed milliseconds from the start time to now
    :param start_dt:
    :return: elapsed time in ms
    """
    return get_diff_millis(start_dt, get_datetime_now())


def get_diff_millis(start_dt, end_dt):
    """
    Calculate millisecond difference between the two datetimes
    :param start_dt:
    :param end_dt:
    :return:
    """
    return convert_timedelta_to_millis(end_dt - start_dt)


def convert_timedelta_to_millis(time_delta):
    """
    Conversion routine

    :param time_delta:
    :return:
    """
    return (time_delta.days * 24 * 60 * 60 + time_delta.seconds) * 1000 + \
           time_delta.microseconds / 1000.0


def on_slurm_cluster():
    """ On the cluster? """
    return 'HPCCLUSTER_INSTITUTE' in os.environ and sys.platform.find('linux') != -1


def on_imac():
    """ On an Mac? """
    return sys.platform.find('darwin') != -1


def get_num_procs():
    """
    Number of CPUs
    :return: number of CPUs
    """

    if on_slurm_cluster():
        return 4
    else:
        import multiprocessing
        return multiprocessing.cpu_count()


def read_yaml(filename):
    """ Read a YAML file """

    with open(filename, 'r') as stream:
        try:
            file_contents = yaml.load(stream)
            return file_contents
        except yaml.YAMLError as exc:
            print(exc)
            exit(1)


def get_file_list(dir_name, suffix_filters=()):
    """
    Get a list of files in a directory

    :param dir_name: name of the directory
    :param suffix_filters: suffices on which to filter
    :return: sorted list of file names
    """

    if not os.path.exists(dir_name):
        raise ValueError("The path '{}' does not exist".format(dir_name))

    if not os.path.isdir(dir_name):
        raise ValueError("The path '{}' is not a directoty".format(dir_name))

    file_names = [f for f in os.listdir(dir_name) if os.path.isfile(os.path.join(dir_name, f))]

    if len(suffix_filters) > 0:
        file_name_set = set()

        for suffix in suffix_filters:
            file_name_set.update({f for f in file_names if f.endswith(suffix)})

        file_names = list(file_name_set)

    return sorted(file_names)


if __name__ == '__main__':
    pass
