""" Optimise the material parameters for a given stomatal geometry """

from __future__ import print_function

import numpy as np
import scipy.optimize as opt
import stomasimulator.febio.febio_sim as fs

from . import stoma_config as sc

c5_factor = 1000.0


def optimise_material_parameters(stoma_cfg):
    """
    Optimise the material parameters for the given stoma to the open shape is attained at the given pressure.

    Parameters
    ----------
    stoma_cfg: sc.StomaConfig
        The stoma configuration

    Returns
    -------
    opt.OptimizeResult:
        The optimisation result
    """

    print('*' * 120)
    print("--> Finding optimum '{}' material parameters for the '{}' stoma...".format(
        stoma_cfg.material_model.label, stoma_cfg.label))
    print('*' * 120)

    optimisation_helper = MaterialParametersOptimisationHelper(stoma_cfg=stoma_cfg)

    # Use full resolution mesh and a different optimiser
    soln = optimisation_helper.do_optimisation()

    return soln


class MaterialParametersOptimisationHelper(object):
    """ Helper class to aid with the optimisation """

    def __init__(self, stoma_cfg, include_lm=False):
        """

        :param stoma_cfg:
        :type stoma_cfg: sc.StomaConfig
        :param include_lm:
        """
        self.stoma_cfg = stoma_cfg
        self.include_lm = include_lm
        self.metrics = None

        # dict of tuple->simulation results (a tuple can be hashed so use it as a key).
        # the optimiser sometimes uses the same parameters so cache the results to avoid unnecessary calculations
        self.cached_results = dict()

    @property
    def material_model(self):
        return self.stoma_cfg.material_model

    def initial_guess(self):
        """
        Get and format the initial guess for the optimisation

        Returns
        -------
        np.ndarray
            The initial guess
        """

        x0 = [self.material_model.isotropic_matrix.c1, self.material_model.isotropic_matrix.c2]

        if not self.material_model.is_isotropic:
            # c5 is scaled in the optimisation function
            x0.append(self.material_model.fibres.c5 / c5_factor)

            if self.include_lm:
                x0.append(self.material_model.fibres.lm)

        return np.asarray(x0)

    def do_optimisation(self):
        """
        Perform the optimisation using SLSQP. Settings tested vs. validation model.

        :return: The optimisation result
        :rtype: OptimizeResult
        """

        print('--> Parameters for optimisation:')
        print('--> Using measurements : {}'.format(self.stoma_cfg.comparison_helper.optimisation_keys))
        print('')

        x0 = self.initial_guess()

        tol, eps = 1e-4, 0.001

        print('--> Using SLSQP with tol={} and eps={}'.format(tol, eps))

        soln = opt.minimize(fun=self.optimise_fn,
                            x0=x0,
                            method='SLSQP',
                            tol=tol,
                            options={'eps': eps})

        print('*' * 120)
        print('--> Optimisation procedure has finished...')
        print(soln)
        print('*' * 120)

        if soln.success:
            print('--> Optimisation succeeded. Result is...')
            self._set_material_parameters(soln.x)
            print('--> {}'.format(self.material_model))
        else:
            print('--> The optimisation failed!')

        print('*' * 120)

        return soln

    def _set_material_parameters(self, x):
        if self.material_model.is_isotropic:
            self.material_model.set_parameters(c1=x[0], c2=x[1])
        else:
            self.material_model.set_parameters(c1=x[0],
                                               c2=x[1],
                                               c5=x[2] * c5_factor,
                                               lm=1.0 if len(x) == 3 else x[3])

        if self.material_model.isotropic_matrix.initial_shear_modulus() <= 0.0:
            print('--> Initial shear modulus is negative! Raw values are {}'.format(x))
            return False

        if not self.material_model.is_isotropic and self.material_model.fibres.c5 < 0.0:
            print('--> Fibre modulus is negative! Raw values are {}'.format(x))
            return False

        return True

    def _bad_metric(self):
        if len(self.cached_results) == 0:
            return 10.0

        return 10 * max(self.cached_results.values())

    def optimise_fn(self, x):
        """
        This method is executed by 'minimize' during the optimisation to evaluate parameter values

        :param x: The current optimisation parameters
        :type x: np.ndarray
        :return: The value that is being minimised
        :rtype: float
        """

        success = self._set_material_parameters(x)
        if not success:
            return self._bad_metric()

        # some iterations are repeated so cache the results to avoid unnecessary iterations
        cached_result_key = tuple(x)
        metric_value = self.cached_results.get(cached_result_key)

        if metric_value is None:
            print('--> Optimiser: {}'.format(self.material_model))

            sim_result = fs.run_simulation(stoma_cfg=self.stoma_cfg,
                                           from_optimiser=True)

            # when the simulation fails we want a non-constant measure for the optimiser to use
            metric_value = sim_result.metric_value if sim_result.success else self._bad_metric()

            self.cached_results[cached_result_key] = metric_value

            print('--> Optimiser: {} - metric={}'.format(self.material_model, metric_value))
        else:
            print('--> Optimiser: {} - metric={} (cached result)'.format(self.material_model, metric_value))

        return metric_value


if __name__ == '__main__':
    raise RuntimeError('Access to this functionality is via the top-level module')
