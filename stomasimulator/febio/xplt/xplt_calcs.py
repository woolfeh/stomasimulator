import stomasimulator.geom.geom_utils as geom


class AttributeCalculator(object):
    """ Abstraction for calculations performed on XPLT state data """

    def __init__(self, prefix, reference_data, dimensionality, lambda_fn=None):
        self.prefix = '' if prefix is None else prefix
        self.reference_data = reference_data
        self.dimensionality = dimensionality
        self.lambda_fn = (lambda x: x) if lambda_fn is None else lambda_fn

    def calculate(self, nid_pt_dict, extras=None):
        """ Perform the calculation

        :param nid_pt_dict: dictionary of an integer 'node id' to a Point object
        :param extras: passed on to the subclass
        :return: a dictionary containing label-result pairs from the calculation
        :rtype: dict
        """

        data = self._calculate(nid_pt_dict, extras)

        if self.dimensionality == 1:
            data = (data,)

        return {k: self.lambda_fn(v) for k, v in zip(self.labels(), data)}

    def _calculate(self, nid_pt_dict, extras):
        """ Calculation implementation - to be overridden in subclasses """
        pass

    def labels(self):
        """ Get the labels for the calculation results """

        suffices = self.calculation_suffices()
        assert len(suffices) == self.dimensionality, 'Error! Data label dimensionality mismatch.'

        fmt_string = '{}{}' if len(self.prefix) == 0 or len(suffices[0]) == 0 else '{}-{}'

        return [fmt_string.format(self.prefix, suffix) for suffix in suffices]

    def calculation_suffices(self):
        """ These suffices are appended to the labels of the calculation result """
        return ['', ] * self.dimensionality


def _get_point(ref_pt, id_pt_dict):
    return id_pt_dict.get(ref_pt) if isinstance(ref_pt, int) else ref_pt


class DistanceCalculator(AttributeCalculator):
    """ Distance between two points """

    def __init__(self, prefix, node_pair, lambda_fn=None):
        node_0 = node_pair[0]
        node_1 = node_pair[1]

        reference_data = (node_0 if node_0.id is None else node_0.id,
                          node_1 if node_1.id is None else node_1.id)

        super(DistanceCalculator, self).__init__(prefix=prefix,
                                                 reference_data=reference_data,
                                                 dimensionality=1,
                                                 lambda_fn=lambda_fn)

    def _calculate(self, nid_pt_dict, extras):
        pt_0 = _get_point(self.reference_data[0], nid_pt_dict)
        pt_1 = _get_point(self.reference_data[1], nid_pt_dict)

        return pt_0.distance(pt_1)


class DirectionalDistanceCalculator(DistanceCalculator):
    """ Signed distance calculator """

    def __init__(self, prefix, node_pair, direction, lambda_fn=None):
        """ Calculate a distance in a specified direction

        :param prefix:
        :param node_pair: two Points - further along 'direction' than node_pair[1] so that 'np[0] - np[1]'
        should be in the direction of 'direction'
        :param direction: the direction vector
        :param lambda_fn:
        """
        super(DirectionalDistanceCalculator, self).__init__(prefix=prefix,
                                                            node_pair=node_pair,
                                                            lambda_fn=lambda_fn)

        self.direction = direction.unit()

    def _calculate(self, nid_pt_dict, extras):
        pt_0 = _get_point(self.reference_data[0], nid_pt_dict)
        pt_1 = _get_point(self.reference_data[1], nid_pt_dict)

        is_in_right_direction = (pt_0 - pt_1) * self.direction > 0.0

        return pt_0.distance(pt_1) if is_in_right_direction else 0.0


class AreaCalculator2D(AttributeCalculator):
    """ Calculate area from a list of points (assumed to be in xy plane) """

    def __init__(self, prefix, boundary_pts, lambda_fn=None):
        super(AreaCalculator2D, self).__init__(prefix=prefix,
                                               reference_data=boundary_pts,
                                               dimensionality=1,
                                               lambda_fn=lambda_fn)

    def _calculate(self, nid_pt_dict, extras):
        updated_pore_pts = [nid_pt_dict[pt.id] for pt in self.reference_data]
        pore_area = geom.calculate_polygon_area(updated_pore_pts)
        return pore_area


class AreaCalculator3D(AttributeCalculator):
    """ Calculate an area from a list of facets """

    def __init__(self, prefix, facet_list):
        super(AreaCalculator3D, self).__init__(prefix=prefix,
                                               reference_data=facet_list,
                                               dimensionality=1)

    def _calculate(self, nid_pt_dict, extras):
        area = geom.calculate_surface_area(nid_pt_dict, self.reference_data)
        return area


class AreaVolumeCalculator(AttributeCalculator):
    """ Perform a combined calculation to get the surface area and volume given a list of facets """

    def __init__(self, prefix, facet_list):
        super(AreaVolumeCalculator, self).__init__(prefix=prefix,
                                                   reference_data=facet_list,
                                                   dimensionality=2)

    def _calculate(self, nid_pt_dict, extras):
        volume, area = geom.calculate_volume_and_area(nid_pt_dict, self.reference_data)
        return area, volume

    def calculation_suffices(self):
        return 'area', 'volume'


class XpltReaderMetrics(object):
    """ Identify the metrics that will be calculated for the XpltReader """

    def __init__(self, comparison_helper=None, is_mesh_calculation_on=False):
        """
        :param comparison_helper: Comparison helper for the stoma
        :type stoma_cfg: sc.ComparisonHelper
        :param is_mesh_calculation_on: Whether to calculate the mesh metrics (or not)
        :type is_mesh_calculation_on: bool
        """
        self.comparison_helper = comparison_helper
        self.is_mesh_calculation_on = is_mesh_calculation_on

    @property
    def is_compare_vs_open_stoma_on(self):
        """
        :return: Whether or not to perform the comparison
        :rtype: bool
        """
        return self.comparison_helper is not None

    def evaluate_metric(self, sim_state):
        """
        Calculate the metric and percent difference vs. each measurement

        :param sim_state: State object holding data from the simulation
        :type sim_state: State
        :return: Each item is a pair comprising a name (key) and its float value
        :rtype: list of tuple
        """

        result = self.comparison_helper.perform_comparison(state_pressure=sim_state.time,
                                                           state_data=sim_state.attributes)

        return result


if __name__ == '__main__':
    pass
