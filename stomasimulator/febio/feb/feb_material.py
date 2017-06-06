"""
The parameters for the nonlinear cell wall materials used by FEBio.

Only the Mooney-Rivlin and the Veronda-Westmann are considered, together with their transversely
isotropic counterparts.
"""

from __future__ import division

from math import (exp, log)

DEFAULT_C5 = 1000.0


class IsotropicMatrix(object):
    """Representation of a two-parameter isotropic matrix"""

    def __init__(self, c1=0.0, c2=0.0):
        self.c1 = c1
        self.c2 = c2

    def __str__(self):
        return '[C1={}, C2={}, K={}]'.format(self.c1, self.c2, self.bulk_modulus())

    def set_parameter_values(self, c1, c2):
        """ Set the parameter values that calibrate the isotropic matrix """
        self.c1 = c1
        self.c2 = c2

    def initial_shear_modulus(self):
        """
        Depends on the material model
        :return:
        """
        pass

    @staticmethod
    def bulk_modulus():
        """
        Use a high bulk modulus, K, to enforce near incompressibility.
        The bulk modulus, K, is set to a constant value of 10,000 MPa for all simulations.

        :return: bulk modulus (a.k.a. penalty parameter in FEBio)
        """

        return 10000.0


class MooneyRivlinMatrix(IsotropicMatrix):
    """ Represent a Mooney-Rivlin matrix """

    def initial_shear_modulus(self):
        """
        G = 2(C1 + C2)
        :return: shear modulus
        """
        return 2.0 * (self.c1 + self.c2)


class VerondaWestmannMatrix(IsotropicMatrix):
    """ Represent a Veronda-Westmann matrix """

    def initial_shear_modulus(self):
        """
        G = C1 * C2
        :return: shear modulus
        """
        return self.c1 * self.c2


class Fibres(object):
    """
    Parameter values for the fibre contribution to the strain energy

    Parameters are: c3, c4, c5, lm (description in FEBio documentation)

    c3 and c4 are set from c5 and lm
    *although*
    they can be set directly by setting them after c5 and lm
    """

    def __init__(self, c3=0.0, c4=0.0, c5=0.0, lm=1.0):
        self.c3 = c3
        self.c4 = c4
        self._c5 = c5
        self._lm = lm
        self._c6 = 0.0

    def __str__(self):
        return '[C5={}, Lm={}]'.format(self.c5, self.lm)

    def set_parameter_values(self, c5=None, lm=1.0):
        """ Set the c5 and lm values """

        self._c5 = c5
        self._lm = lm

        self._update()

    @property
    def c5(self):
        """ Get the fibre modulus """
        return self._c5

    @property
    def c6(self):
        """ Get the c6 parameter """
        return self._c6

    @property
    def lm(self):
        """ Return the lambda_m property """
        return self._lm

    def _update(self):
        self._set_implied_c3_c4()
        self._set_c6()

    def _set_implied_c3_c4(self):
        """
        To reduce parameter space we impose cty on the derivative of the fibre strain energy to get
        c3 and c4 from c5 and lm

        Since c3 * c4 is the slope of the stress curve at lambda = 1+, we set (for convenience)
        c3 * c4 = c5 / 100 (i.e. 1% of c5)

        Slope is cts at lambda = lm if:
          c3 * c4 * exp( c4*(lm-1) ) = c5
        so
          c4 = ln( c5/(c3*c4) ) / (lm-1) = ln( 100 ) / (lm - 1)

        This leaves c5 and lm as the free fibre parameters - a two parameter space for the fibres

        Note: if lm == 1 then there is no exponential stress region so set nominal values
        """

        if self.c5 is None or self.lm is None:
            c3, c4 = None, None
        elif self.lm == 1.0:
            # set nominal values
            c3, c4 = 0.0, 0.0
        else:
            divisor = 100.0
            c3c4 = self.c5 / divisor

            c4 = 1.0 if self.lm == 1.0 else log(divisor) / (self.lm - 1.0)
            c3 = c3c4 / c4

        self.c3, self.c4 = c3, c4

    def _set_c6(self):
        self._c6 = self.fibre_strain_energy(self.lm) - self.c5 * self.lm

    def fibre_strain_energy(self, l_stretch):
        """ Calculate the strain energy of the fibres
        :param l_stretch:
        :return:
        """
        if l_stretch <= 1.0:
            # compressed region - no energy
            return 0.0

        # Note: this range should be '< lm' according to FEBio but we use '<=' to
        #       make setting c6 easier -> there's no difference because it's cts.
        if l_stretch <= self.lm:
            # exponential energy
            return self.c3 * (exp(self.c4 * (l_stretch - 1.0)) - 1.0)

        # linear energy
        return self.c5 * l_stretch + self.c6


class MaterialModel(object):
    """ Representation of a material model """

    def __init__(self, isotropic_matrix, fibres=None):
        """
        Material model for FEBio - aggregates matrix and fibre objects

        :type isotropic_matrix: IsotropicMatrix
        :type fibres: Fibres

        :param isotropic_matrix:
        :param fibres:
        :return:
        """
        self.isotropic_matrix = isotropic_matrix
        self.fibres = fibres

    def __str__(self):
        if self.is_isotropic:
            return "['{}', matrix:{}]".format(self.febio_name, self.isotropic_matrix)
        else:
            return "['{}', matrix:{}, fibres:{}]".format(self.febio_name,
                                                         self.isotropic_matrix,
                                                         self.fibres)

    @property
    def febio_name(self):
        """ The name in FEBio """
        return ''

    @property
    def label(self):
        """ Label for the material model """
        return ''

    @property
    def is_isotropic(self):
        """ Is the material isotropic?
        :return: True when fibres are absent
        """
        return self.fibres is None

    def set_parameters(self, c1, c2, c5=None, lm=1.0):
        """ Set material parameters """

        self.isotropic_matrix.set_parameter_values(c1=c1, c2=c2)

        if not self.is_isotropic:
            self.fibres.set_parameter_values(c5=c5 if c5 else DEFAULT_C5,
                                             lm=lm)


class MooneyRivlinMaterialModel(MaterialModel):
    """ Represent a Mooney-Rivlin material """

    def __init__(self, c1=5.0, c2=5.0):
        MaterialModel.__init__(self, MooneyRivlinMatrix(c1=c1, c2=c2))

    @property
    def febio_name(self):
        return 'Mooney-Rivlin'

    @property
    def label(self):
        return 'MR'


class TransIsoMooneyRivlinMaterialModel(MaterialModel):
    """ Represent a transversely-isotropic Mooney-Rivlin material """

    def __init__(self, c1=5.0, c2=5.0, c5=DEFAULT_C5, lm=1.0):
        MaterialModel.__init__(self, MooneyRivlinMatrix(c1, c2), Fibres(c5=c5, lm=lm))

    @property
    def febio_name(self):
        return 'trans iso Mooney-Rivlin'

    @property
    def label(self):
        return 'tiMR'


class VerondaWestmannMaterialModel(MaterialModel):
    """ Represent a Veronda-Westmann material """

    def __init__(self, c1=5.0, c2=5.0):
        MaterialModel.__init__(self, VerondaWestmannMatrix(c1=c1, c2=c2))

    @property
    def febio_name(self):
        return 'Veronda-Westmann'

    @property
    def label(self):
        return 'VW'


class TransIsoVerondaWestmannMaterialModel(MaterialModel):
    """ Represent a transversely-isotropic Veronda-Westmann material """

    def __init__(self, c1=5.0, c2=5.0, c5=DEFAULT_C5, lm=1.0):
        MaterialModel.__init__(self, VerondaWestmannMatrix(c1=c1, c2=c2), Fibres(c5=c5, lm=lm))

    @property
    def febio_name(self):
        return 'trans iso Veronda-Westmann'

    @property
    def label(self):
        return 'tiVW'


MATERIAL_MODEL_CLASSES = dict()

MATERIAL_MODEL_LABELS = [MooneyRivlinMaterialModel().label,
                         VerondaWestmannMaterialModel().label,
                         TransIsoMooneyRivlinMaterialModel().label,
                         TransIsoVerondaWestmannMaterialModel().label]


def _populate_material_model_dict():
    if len(MATERIAL_MODEL_CLASSES) > 0:
        return

    for mm in (MooneyRivlinMaterialModel,
               VerondaWestmannMaterialModel,
               TransIsoMooneyRivlinMaterialModel,
               TransIsoVerondaWestmannMaterialModel):
        MATERIAL_MODEL_CLASSES[mm().label] = mm


def lookup_material_model(material_model_label):
    """ Link a label to a new instance of the material model

    :param material_model_label:
    :return: MaterialModel object
    :rtype: MaterialModel
    """
    _populate_material_model_dict()
    mm_cls = MATERIAL_MODEL_CLASSES[material_model_label]
    return mm_cls()


if __name__ == '__main__':
    raise NotImplementedError('Access to this functionality is via the top-level module')
