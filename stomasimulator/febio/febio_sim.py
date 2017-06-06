""" Run FEBio simulation """

from __future__ import print_function

import os
import subprocess

import stomasimulator.config as config
import stomasimulator.stomata.stoma_config as sc
import stomasimulator.utils.utils as ut

from .feb import feb_writer as fw
from .xplt import (xplt_calcs as xcalc,
                   xplt_reader as xr)


def generate_filename_base(stoma_cfg):
    """
    :param stoma_cfg:
    :type stoma_cfg: sc.StomaConfig
    :param filename:
    :return: basename for the FEBio files
    :rtype: str
    """

    # Build the file name using the stoma config label, some geometry information and
    # then the material parameters

    fn_parts = [stoma_cfg.label.strip(),
                '{}gc'.format(stoma_cfg.num_gcs),
                'd_{}'.format(stoma_cfg.wall_thickness.dorsal),
                'v_{}'.format(stoma_cfg.wall_thickness.ventral)]

    if stoma_cfg.aspect_ratio != 1.0:
        fn_parts.append('ar_{}'.format(stoma_cfg.aspect_ratio))

    mm = stoma_cfg.material_model

    fn_parts.append(mm.label)
    fn_parts.append('c1_{:.4f}'.format(mm.isotropic_matrix.c1))
    fn_parts.append('c2_{:.4f}'.format(mm.isotropic_matrix.c2))

    if not mm.is_isotropic:
        fn_parts.append('c5_{:.2f}'.format(mm.fibres.c5))
        fn_parts.append('lm_{:.2f}'.format(mm.fibres.lm))

    filename_base = ','.join(fn_parts)

    return filename_base


class FEBioSimulationResult(object):
    """ Store file names and results for the FEBio simulation """

    def __init__(self, stoma_cfg, filename_base=None, from_optimiser=False):
        """
        Initialise the simulation result
        :param stoma_cfg: Stoma configuration object
        :param filename_base: base filename - no extension
        :param from_optimiser: whether called from the optimiser
        """

        if filename_base:
            self.filename_base = filename_base
        else:
            self.filename_base = generate_filename_base(stoma_cfg)

        self.from_optimiser = from_optimiser

        self.final_state = None
        self.metric_value = None
        self.all_diffs = None

    def __str__(self):
        return '[feb file: {}, xplt file: {}, metric: {}, final state: {}'.format(
            self.feb_filename, self.xplt_filename, self.metric_value, self.final_state)

    @property
    def feb_filename(self):
        """ Return the name of the feb file """
        return '{}.feb'.format('temp_optimiser' if self.from_optimiser else self.filename_base)

    @property
    def xplt_filename(self):
        """ Return the name of the xplt file """
        return '{}.xplt'.format('temp_optimiser' if self.from_optimiser else self.filename_base)

    @property
    def results_filename(self):
        """ Return the name of the xplt file """
        return '{}.stats.txt'.format(self.filename_base)

    @property
    def success(self):
        """
        Was the simulation successful?
        :return: True when the final state and the metric value are both set
        """
        return self.final_state is not None and self.metric_value is not None


def run_febio_process(feb_filename, xplt_filename):
    """ Run the FEBio simulation """

    command = [config.APP_CONFIG.febio_exe,
               '-silent',
               '-i', feb_filename,
               '-p', xplt_filename]

    print('--> Running simulation with command: {}'.format(' '.join(command)))

    # set number of OpenMP threads - to make sure it is set
    my_env = os.environ.copy()
    my_env['OMP_NUM_THREADS'] = str(ut.get_num_procs())

    # run the FEBio process - takes several minutes...
    subprocess.call(command, env=my_env)


def process_xplt_file(sim_result, metrics):
    """
    Process the XPLT file produced by the FEBio simulation

    :param sim_result:
    :type sim_result: FEBioSimulationResult
    :param metrics:
    :type metrics: xcls.XpltReaderMetrics
    :return:
    :rtype: FEBioSimulationResult
    """

    # process the file and get the last state
    sim_result.final_state = xr.process_xplt_file(xplt_filename=sim_result.xplt_filename,
                                                  results_filename=sim_result.results_filename,
                                                  metrics=metrics)

    if sim_result.final_state is not None:
        sim_result.metric_value = sim_result.final_state.attributes['metric']

        # For the LM optimiser
        sim_result.all_diffs = sim_result.final_state.attributes['all-diffs']

    return sim_result


def run_simulation(stoma_cfg,
                   do_mesh_calculation=False,
                   from_optimiser=False,
                   filename_base=None,
                   only_generate_feb_file=False):
    """
    :param stoma_cfg:
    :type stoma_cfg: sc.StomaConfig
    :param do_mesh_calculation:
    :param from_optimiser:
    :param filename_base:
    :param only_generate_feb_file:
    :return:
    :rtype: FEBioSimulationResult
    """

    # set up the simulation result
    sim_result = FEBioSimulationResult(stoma_cfg=stoma_cfg,
                                       from_optimiser=from_optimiser,
                                       filename_base=filename_base)

    fw.output_febio_xml(stoma_cfg=stoma_cfg,
                        filename=sim_result.feb_filename,
                        perform_mesh_checks=do_mesh_calculation)

    if only_generate_feb_file:
        return sim_result

    run_febio_process(feb_filename=sim_result.feb_filename,
                      xplt_filename=sim_result.xplt_filename)

    metrics = xcalc.XpltReaderMetrics(comparison_helper=stoma_cfg.comparison_helper,
                                      is_mesh_calculation_on=do_mesh_calculation)

    process_xplt_file(sim_result, metrics)

    return sim_result


if __name__ == '__main__':
    raise NotImplementedError('Access to this functionality is via the top-level module')
