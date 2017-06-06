#!/usr/bin/env python

"""
Brings all the different modules together so that they can be run from one place.
The argument parsers are in each dispatched method.
"""

from __future__ import (absolute_import, print_function)

import argparse
import sys
from math import sqrt

from sortedcontainers import SortedDict

import stomasimulator.febio.feb.feb_material as fm
import stomasimulator.febio.febio_sim as febsim
import stomasimulator.febio.xplt.xplt_calcs as xcalc
import stomasimulator.febio.xplt.xplt_reader as xr
import stomasimulator.stomata.opt_material as mc
import stomasimulator.stomata.stoma_config as sc

from . import config

# keys must match the names of class methods that perform that action
COMMAND_INFO = SortedDict(cfg="Output the configuration settings and model labels used by stomasimulator",
                          feb='Output a FEB file for the specified stoma model (the simulation can be run)',
                          xplt='Process an XPLT file',
                          opt='Optimise the material parameters for a given stoma model', )


class Dispatcher(object):
    """
    Read arguments and dispatch the appropriate function
    """

    def __init__(self):
        """ Initialise object """

        parser = argparse.ArgumentParser(description='Execute functions within the stomata biomechanics suite',
                                         usage='\n'.join(self.usage_string()))

        parser.add_argument('command', help='Subcommand to run')

        # use argv to get the command (use a range so that an exception is not raised if there are no arguments)
        args = parser.parse_args(sys.argv[1:2])

        if not hasattr(self, args.command):
            print('Unrecognised command')
            parser.print_help()
            exit(1)

        # use dispatch pattern to invoke method with same name
        description = COMMAND_INFO[args.command]
        prog = 'stomasimulator {}'.format(args.command)

        command_parser = argparse.ArgumentParser(prog=prog, description=description)

        getattr(self, args.command)(command_parser, sys.argv[2:])

    @staticmethod
    def usage_string():
        """
        Generate the 'how to use me' string
        :return:
        """

        usage_str = ['stomasimulator <command> [<args>]', '', 'The commands are:']
        for key, value in COMMAND_INFO.items():
            usage_str.append('  {:12s} {}'.format(key, value))
        usage_str.append('')

        return usage_str

    @staticmethod
    def _generate_stoma_config(model_label,
                               whole_stoma,
                               material_model_label=None,
                               verbose=True,
                               optimise_vs_open=True):
        """
        Create a StomaConfig object given the label

        :param model_label:
        :param whole_stoma:
        :return: stoma configuration object
        :rtype: sc.StomaConfig
        """

        if verbose:
            print("--> Looking for a '{}' stoma...".format(model_label), end='')

        try:
            yaml_file = config.get_file_contents_for_model_label(model_label)
        except KeyError:
            print('{} stoma not found'.format(model_label))
            exit(1)

        if verbose:
            print('found')

        scfg = sc.generate_stoma_model(yaml_file=yaml_file,
                                       model_label=model_label,
                                       material_label=material_model_label,
                                       whole_stoma=whole_stoma,
                                       optimise_vs_open=optimise_vs_open)

        return scfg

    ### The dispatch methods...

    @classmethod
    def cfg(cls, parser, args):
        """
        Output the configuration settings

        :param parser:
        :type: argparse.ArgumentParser
        :param args: ignored
        :return:
        """

        parser.add_argument('-f', '--model-filter', help='Search for models with names containing the given string')

        args = parser.parse_args(args)

        config.print_model_labels(args.model_filter)

    @classmethod
    def feb(cls, parser, args):
        """
        Generate the XML configuration file for an FEBio simulation

        :param parser:
        :type: argparse.ArgumentParser
        :param args:
        :return:
        """

        # mandatory
        parser.add_argument('model_label', help='Label of a model stoma in the config. file')
        parser.add_argument('material', help='Label for the material model', choices=fm.MATERIAL_MODEL_LABELS)

        # optional
        parser.add_argument('-w', '--whole-stoma', action='store_true',
                            help='If set then output the whole stoma (default is one guard cell)')
        parser.add_argument('-s', '--simulation', action='store_true',
                            help='Run simulation for the model')

        # for parameter adjustment
        parser.add_argument('--matrix-factor', nargs='?', type=float, help='Multiply G0 by this factor')
        parser.add_argument('--fibre-factor', nargs='?', type=float, help='Multiply C5 by this factor')

        args = parser.parse_args(args)

        # args.simulation = False

        # Set the stoma configuration
        scfg = cls._generate_stoma_config(model_label=args.model_label,
                                          material_model_label=args.material,
                                          whole_stoma=args.whole_stoma)

        if args.matrix_factor is not None or args.fibre_factor is not None:
            mm = scfg.material_model

            filename_parts = [args.model_label, args.material, 'sens']

            # TODO Move parameter adjustment to the 'febmat' module

            if args.matrix_factor is not None:
                print('--> Using matrix factor', args.matrix_factor)

                filename_parts.append('G0x{}'.format(args.matrix_factor))

                # In the Mooney-Rivlin model, G_0 = 2 * ( C1 + C2 ), so to test the sensitivity of G to a change,
                # lambda,multiply by lambda:
                #   lambda * G_0 == lambda * 2 * ( C1 + C2 ) = 2 * ( lambda * C1 + lambda * C2 )
                #
                # In the Veronda-Westmann model, G_0 = C1 * C2, so to test the sensitivity of G to a change, lambda,
                # multiply by the sqrt of lambda:
                #   lambda * G == lambda * C1 * C2 = ( sqrt( lambda ) * C1 ) * ( sqrt( lambda ) * C2 )

                c_fac = args.matrix_factor
                if mm.label == fm.TransIsoVerondaWestmannMaterialModel().label:
                    c_fac = sqrt(args.matrix_factor)

                scfg.material_model.isotropic_matrix.set_parameter_values(c1=mm.isotropic_matrix.c1 * c_fac,
                                                                          c2=mm.isotropic_matrix.c2 * c_fac)

            if args.fibre_factor is not None:
                print('--> Using fibre factor', args.fibre_factor)

                filename_parts.append('C5x{}'.format(args.fibre_factor))

                scfg.material_model.fibres.set_parameter_values(c5=mm.fibres.c5 * args.fibre_factor,
                                                                lm=mm.fibres.lm)

            filename_base = ','.join(filename_parts)
        else:
            filename_base = None

        if args.simulation:
            febsim.run_simulation(stoma_cfg=scfg, filename_base=filename_base)
        else:
            # just generate the FEB file
            febsim.run_simulation(stoma_cfg=scfg, filename_base=filename_base, only_generate_feb_file=True)

    @classmethod
    def opt(cls, parser, args):
        """
        For the given stomatal configuration optimise the material parameters

        :param parser:
        :param args:
        :return:
        """

        # mandatory
        parser.add_argument('model_label', help='Label of a model stoma')
        parser.add_argument('material', help='Label for the material model', choices=fm.MATERIAL_MODEL_LABELS)

        # optional
        parser.add_argument('-w', '--whole-stoma', action='store_true',
                            help='If set then use the whole stoma (default is one guard cell)')

        args = parser.parse_args(args)

        # Use the stoma configuration
        scfg = cls._generate_stoma_config(model_label=args.model_label,
                                          material_model_label=args.material,
                                          whole_stoma=args.whole_stoma)

        mc.optimise_material_parameters(stoma_cfg=scfg)

    @classmethod
    def xplt(cls, parser, args):
        """
        Process an XPLT file from FEBio
        :param parser:
        :param args:
        :return:
        """

        parser.add_argument('xplt_file_name', help='Name of the xplt file to process')

        parser.add_argument('-l', '--model-label', help='Label of a model stoma in the config. file')
        parser.add_argument('-m', '--mesh-metrics', help='Use VTK to calculate the mesh metrics (at start and end)',
                            action='store_true')
        parser.add_argument('-r', action='store_true', help='Output the raw conversion of the xplt file and exit')
        parser.add_argument('-v', action='store_true', help='Verbose flag')

        args = parser.parse_args(args)

        stoma_cfg = None
        if args.model_label is not None:
            # create a StomaConfig object

            whole_stoma = '2gc' in args.xplt_file_name

            stoma_cfg = cls._generate_stoma_config(model_label=args.model_label,
                                                   whole_stoma=whole_stoma)

        metrics = xcalc.XpltReaderMetrics(comparison_helper=None if stoma_cfg is None else stoma_cfg.comparison_helper,
                                          is_mesh_calculation_on=args.mesh_metrics)

        xr.process_xplt_file(xplt_filename=args.xplt_file_name,
                             verbose=args.v,
                             output_raw_data=args.r,
                             metrics=metrics)

        return


if __name__ == '__main__':
    pass
