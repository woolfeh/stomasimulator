""" Represent planar and 3D ellipses """

from __future__ import print_function

from math import pi

import numpy as np

import scipy as sp
import scipy.optimize as opt
import scipy.special as spsp

from . import point as p
from . import geom_utils


def calculate_cost_sint(t):
    """
    Calculate cos(t) and sin(t) for the given value of t
    :param t:
    :return: 2-tuple of ( cos(t), sin(t) )
    """
    if abs(2 * t - pi) < 1e-10:
        cost, sint = 0.0, 1.0
    else:
        cost, sint = np.cos(t), np.sin(t)

    return cost, sint


class Ellipse(object):
    """
    Represent a planar (in x-y plane) ellipse
    """

    def __init__(self, semi_x_axis, semi_y_axis):
        """ """
        self.__semi_x_axis = float(semi_x_axis)
        self.__semi_y_axis = float(semi_y_axis)

        self.__m = None
        self.__circumference = None

    def __str__(self):
        """ """
        return '[Ellipse: a={}, b={}]'.format(self.semi_x_axis, self.semi_y_axis)

    @property
    def semi_x_axis(self):
        """ """
        return self.__semi_x_axis

    @property
    def semi_y_axis(self):
        """ """
        return self.__semi_y_axis

    @property
    def m(self):
        """
        Returns the value, m, required by the SciPy elliptic integral functions.
        The definition is
            m = 1 - (minor axis/major axis)^2
        and here we take major to be x
        The ellipse eccentricity, e, and m are related by m = e^2
        """
        if self.__m is None:
            self.__m = 1.0 - (self.semi_y_axis / self.semi_x_axis) ** 2

        return self.__m

    @property
    def circumference(self):
        """
        Total perimeter of an ellipse is, L = 4 a E(m), where a is the major axis and E(m) is the
        complete elliptical integral

        :return:
        """
        if self.__circumference is None:
            self.__circumference = 4 * self.semi_x_axis * spsp.ellipe(self.m)

        return self.__circumference

    def calculate_ellipse_xy_for_t(self, t):
        """
        Ellipse is (x/a)**2 + (y/b)**2 = 1 which can be parameterised using (x,y)=(a cos(t), b sin(t))
        Note that 't' is the parametric angle and is *not* the polar angle, 'p' [ tan(p) = y/x = (b/a) tan(t) ]

        :param t: the parametric angle
        :return: x, y coordinates
        """

        cost, sint = calculate_cost_sint(t)

        x, y = self.semi_x_axis * cost, self.semi_y_axis * sint

        return x, y

    def calculate_ellipse_pt_for_t(self, t):
        """
        Returns a Point object (in the xy plane) for the parametric angle theta

        :param t: parametric angle
        :return:
        """
        x, y = self.calculate_ellipse_xy_for_t(t)

        return p.Point(x, y, 0.0)

    def calculate_ellipse_t_for_pt(self, pt):
        """
        Returns the parametric angle for the point

        :param pt:
        :return:
        """

        from math import acos, asin

        if abs(pt.x) < abs(pt.y):
            t = acos(pt.x / self.semi_x_axis)
        else:
            t = asin(pt.y / self.semi_y_axis)

        return t

    def calculate_ellipse_normal_for_t(self, t, outward=True):
        """
        Point on the ellipse is ( a cos(t), b sin(t) ) so the outward normal is ( b cos(t), a sin(t) )

        :param t: parametric angle
        :param outward: whether to make it outward pointing or inward
        :return: unit normal at t
        """

        cost, sint = calculate_cost_sint(t)

        norm = p.Point(self.semi_y_axis * cost, self.semi_x_axis * sint, 0.0).unit()

        return norm if outward else -norm

    def calculate_ellipse_x(self, y):
        """
        Calculate the x coordinate (in quadrant 1) for a given y value
        :param y:
        :return: x
        """

        assert y is not None and y >= 0.0, 'Invalid value for y (y={})'.format(y)

        x = self.semi_x_axis * sp.sqrt(1.0 - (y / self.semi_y_axis) ** 2)
        return x

    def calculate_ellipse_y(self, x):
        """
        Calculate the y coordinate (in quadrant 1) for a given x value
        :param x:
        :return: y
        """

        assert x is not None and x >= 0.0, 'Invalid value for x (x={})'.format(x)

        y = self.semi_y_axis * sp.sqrt(1.0 - (x / self.semi_x_axis) ** 2)

        return y

    def calculate_curvature_for_t(self, t):
        """

        :param t:
        :return:
        """

        cost, sint = calculate_cost_sint(t)

        numer = self.semi_x_axis * self.semi_y_axis
        denom = ((self.semi_x_axis * sint) ** 2 + (self.semi_y_axis * cost) ** 2) ** 1.5

        k = numer / denom

        return k


class SuperEllipse(object):
    """
    Represent a planar (in x-y plane) ellipse
    """

    def __init__(self, semi_x_axis, semi_y_axis, r):
        """ """
        self.semi_x_axis = float(semi_x_axis)
        self.semi_y_axis = float(semi_y_axis)
        self.r = r

        self.__m = None
        self.__circumference = None

    def __str__(self):
        """ """
        return '[Ellipse: a={}, b={}, r={}]'.format(self.semi_x_axis, self.semi_y_axis, self.r)

    def calculate_xy_for_t(self, t):
        """
        Ellipse is (x/a)**r + (y/b)**r = 1 which can be parameterised using
            (x,y) = ( a cos(t)**(2/r), b sin(t)**(2/r))
        Note that 't' is the parametric angle and is *not* the polar angle, 'p' [ tan(p) = y/x = (b/a) tan(t) ]

        :param t: the parametric angle
        :return: x, y coordinates
        """

        cost, sint = calculate_cost_sint(t)

        x, y = self.semi_x_axis * cost ** (2.0 / self.r), self.semi_y_axis * sint ** (2.0 / self.r)

        return x, y

    def calculate_ellipse_pt_for_t(self, t):
        """
        Returns a Point object (in the xy plane) for the parametric angle theta

        :param t: parametric angle
        :return:
        """
        x, y = self.calculate_xy_for_t(t)

        return p.Point(x, y, 0.0)

    def calculate_t_for_pt(self, pt):
        """
        Returns the parametric angle for the point

        :param pt:
        :return:
        """

        from math import acos, asin

        if abs(pt.x) < abs(pt.y):
            t = acos((pt.x / self.semi_x_axis) ** (self.r / 2))
        else:
            t = asin((pt.y / self.semi_y_axis) ** (self.r / 2))

        return t


class Ellipse3D(object):
    """
    Any 3D ellipse can be rotated and translated to a planar x-y ellipse, so this class wraps the above
    Ellipse class with calculations performed there and the necessary rotation and translation performed here.
    """

    def __init__(self, axis_pt_1, axis_pt_2, centre_pt):
        """
        An 3D ellipse

        :param axis_pt_1: Point object - a point at one axis
        :param axis_pt_2: Point object - a point at the other axis
        :param centre_pt: Point object - Centre of the ellipse, which doubles as the translation vector
        :return:
        """

        pt_1 = axis_pt_1 - centre_pt
        pt_2 = axis_pt_2 - centre_pt

        # make sure the centre point is central - scalar product of vectors should be zero
        assert abs(pt_1 * pt_2) < 1e-10, 'The given centre point is not central'

        self._axis_pt_1 = axis_pt_1
        self._axis_pt_2 = axis_pt_2
        self._centre_point = centre_pt

        # translate the points - make pt_1 the origin
        pt_1 = pt_1.xyz
        pt_2 = pt_2.xyz

        # put point 1 on the x-z plane by rotating about the z-axis
        angle = np.arctan2(pt_1[1], pt_1[0])
        rot_mat_z = geom_utils.rotation_matrix([0, 0, 1], -angle)
        pt_1, pt_2 = np.dot(rot_mat_z, pt_1), np.dot(rot_mat_z, pt_2)

        # rotate about the y axis to put axis point 1 on the x axis
        angle = np.arctan2(pt_1[0], pt_1[2])
        rot_mat_y = geom_utils.rotation_matrix([0, 1, 0], pi / 2 - angle)
        pt_1, pt_2 = np.dot(rot_mat_y, pt_1), np.dot(rot_mat_y, pt_2)

        # rotate about the x axis to put axis point 2 on the y axis (- want to send it clockwise)
        angle = np.arctan2(pt_2[2], pt_2[1])
        rot_mat_x = geom_utils.rotation_matrix([1, 0, 0], -angle)
        pt_1, pt_2 = np.dot(rot_mat_x, pt_1), np.dot(rot_mat_x, pt_2)

        # assert pt_1[0] > 0.0 and pt_2[1] > 0.0 and pt_1[2] == 0.0 and pt_2[2] == 0.0

        self._ellipse = Ellipse(pt_1[0], pt_2[1])

        self._rotation_matrix = np.dot(rot_mat_x, np.dot(rot_mat_y, rot_mat_z))
        self._inverse_rotation = np.transpose(self._rotation_matrix)

    @property
    def axis_pt_1(self):
        """ """
        return self._axis_pt_1

    @property
    def axis_pt_2(self):
        """ """
        return self._axis_pt_2

    @property
    def centre_pt(self):
        """ """
        return self._centre_point

    def _convert_local_pt_to_global_pt(self, local_pt):
        """
        Rotate and translate the local point to the global frame of reference
        :param local_pt:
        :return:
        """
        rot_xyz = np.dot(self._inverse_rotation, local_pt.xyz)
        global_pt = self.centre_pt + p.from_triad(rot_xyz)
        return global_pt

    def _convert_global_pt_local_pt(self, global_pt):
        """
        Rotate and translate the global point to the local frame of reference
        :param global_pt:
        :return:
        """
        tran_pt = global_pt - self.centre_pt
        local_pt = np.dot(self._rotation_matrix, tran_pt.xyz)
        return p.from_triad(local_pt)

    def calculate_equi_arc_spaced_t(self, num):
        """
        Calculate the t-values (parametric angles) for the quarter ellipse.
        :param num:
        :return:
        """
        t_values = calculate_equi_arc_spaced_t(self._ellipse, num)
        return t_values

    def calculate_ellipse_pt_for_t(self, t):
        """

        :param t:
        :return:
        """
        local_pt = self._ellipse.calculate_ellipse_pt_for_t(t)
        return self._convert_local_pt_to_global_pt(local_pt)

    def calculate_equi_arc_spaced_pts(self, num):
        """

        :param num:
        :return:
        """
        t_values = self.calculate_equi_arc_spaced_t(num)
        pts = [self.calculate_ellipse_pt_for_t(t) for t in t_values]
        return pts

    def calculate_equi_t_spaced_pts(self, num):
        """
        Calculate a list of Points for the equi-spaced t-values
        :param num:
        :return:
        """
        local_pts = calculate_pts_for_equi_spaced_t(self._ellipse, num)
        pts = [self._convert_local_pt_to_global_pt(lp) for lp in local_pts]
        return pts


class SuperEllipse3D(object):
    """
    Any 3D ellipse can be rotated and translated to a planar x-y ellipse, so this class wraps the above
    Ellipse class with calculations performed there and the necessary rotation and translation performed here.
    """

    def __init__(self, axis_pt_1, axis_pt_2, centre_pt, r):
        """
        An 3D ellipse

        :param axis_pt_1: Point object - a point at one axis
        :param axis_pt_2: Point object - a point at the other axis
        :param centre_pt: Point object - Centre of the ellipse, which doubles as the translation vector
        :return:
        """

        pt_1 = axis_pt_1 - centre_pt
        pt_2 = axis_pt_2 - centre_pt

        # make sure the centre point is central - scalar product of vectors should be zero
        assert abs(pt_1 * pt_2) < 1e-10, 'The given centre point is not central'

        self.axis_pt_1 = axis_pt_1
        self.axis_pt_2 = axis_pt_2
        self.centre_pt = centre_pt
        self.r = r

        # translate the points - make pt_1 the origin
        pt_1 = pt_1.xyz
        pt_2 = pt_2.xyz

        # put point 1 on the x-z plane by rotating about the z-axis
        angle = np.arctan2(pt_1[1], pt_1[0])
        rot_mat_z = geom_utils.rotation_matrix([0, 0, 1], -angle)
        pt_1, pt_2 = np.dot(rot_mat_z, pt_1), np.dot(rot_mat_z, pt_2)

        # rotate about the y axis to put axis point 1 on the x axis
        angle = np.arctan2(pt_1[0], pt_1[2])
        rot_mat_y = geom_utils.rotation_matrix([0, 1, 0], pi / 2 - angle)
        pt_1, pt_2 = np.dot(rot_mat_y, pt_1), np.dot(rot_mat_y, pt_2)

        # rotate about the x axis to put axis point 2 on the y axis (- want to send it clockwise)
        angle = np.arctan2(pt_2[2], pt_2[1])
        rot_mat_x = geom_utils.rotation_matrix([1, 0, 0], -angle)
        pt_1, pt_2 = np.dot(rot_mat_x, pt_1), np.dot(rot_mat_x, pt_2)

        # assert pt_1[0] > 0.0 and pt_2[1] > 0.0 and pt_1[2] == 0.0 and pt_2[2] == 0.0

        self._ellipse = SuperEllipse(pt_1[0], pt_2[1], r)

        self._rotation_matrix = np.dot(rot_mat_x, np.dot(rot_mat_y, rot_mat_z))
        self._inverse_rotation = np.transpose(self._rotation_matrix)

    def _convert_local_pt_to_global_pt(self, local_pt):
        """
        Rotate and translate the local point to the global frame of reference
        :param local_pt:
        :return:
        """
        rot_xyz = np.dot(self._inverse_rotation, local_pt.xyz)
        global_pt = self.centre_pt + p.from_triad(rot_xyz)
        return global_pt

    def _convert_global_pt_local_pt(self, global_pt):
        """
        Rotate and translate the global point to the local frame of reference
        :param global_pt:
        :return:
        """
        tran_pt = global_pt - self.centre_pt
        local_pt = np.dot(self._rotation_matrix, tran_pt.xyz)
        return p.from_triad(local_pt)

    def calculate_ellipse_pt_for_t(self, t):
        """

        :param t:
        :return:
        """
        local_pt = self.super_ellipse.calculate_pt_for_t(t)
        return self._convert_local_pt_to_global_pt(local_pt)

    def calculate_equi_t_spaced_pts(self, num):
        """
        Calculate a list of Points for the equi-spaced t-values
        :param num:
        :return:
        """
        local_pts = calculate_pts_for_equi_spaced_t(self._ellipse, num)
        pts = [self._convert_local_pt_to_global_pt(lp) for lp in local_pts]
        return pts


class EllipticIntegralHelper(object):
    """
    Help
    """

    @classmethod
    def convert_phi_to_t(cls, phi_or_t):
        """
        Calculate the parametric angle, t, from the elliptic integral amplitude, phi (or vice versa)
        Note: both phi and t are parametric angles - NOT polar angles

        :param phi_or_t: amplitude (phi is measured from the y-axis)
        :return: theta or phi
        """
        th_ph = pi / 2 - phi_or_t
        return 0.0 if abs(th_ph) < 1e-10 else th_ph

    @classmethod
    def find_phi_for_arc_length(cls, ellipse, arc):
        """
        Computes the elliptic integral angle, phi, for the given arc-length - measured from the y-axis
        Only considers the ellipse quarter that lies in the first quadrant

        :param ellipse:
        :param arc:
        :return:
        """

        if arc == 0.0:
            phi = 0.0
        else:
            def opt_fn(x):
                """ """
                return ellipse.semi_x_axis * spsp.ellipeinc(x, ellipse.m) - arc

            # optimiser sometimes needs some margin so add 0.01 at each end
            phi_tuple = opt.brentq(opt_fn, -0.01, pi / 2 + 0.01, disp=True)

            phi = phi_tuple[0] if isinstance(phi_tuple, tuple) else phi_tuple

        return phi


def calculate_ellipse_arc(ellipse, t):
    """

    :param ellipse:
    :param t:
    :return:
    """

    phi = EllipticIntegralHelper.convert_phi_to_t(t)

    ellipse_arc_from_y_axis = ellipse.semi_x_axis * spsp.ellipeinc(phi, ellipse.m)

    ellipse_arc = ellipse.circumference / 4 - ellipse_arc_from_y_axis

    return ellipse_arc


def calculate_equi_arc_spaced_t(ellipse, num, from_t=None):
    """
    Calculate values of t that divide quarter of the ellipse into arcs of equal length
    Note: The arc length variables in this function refer to the arc length along the ellipse quarter
          from the y-axis

    :param ellipse:
    :param num:
    :param from_t:
    :return:
    """

    def convert_t_to_theta(_t):
        """
        tan( theta ) = b sin t / a cost t = (b/a) * tan t

        :param _t:
        :return:
        """

        x, y = ellipse.calculate_ellipse_xy_for_t(_t)

        return np.arctan2(y, x)

    def convert_theta_to_t(_theta):
        """

        :param _theta:
        :return:
        """
        return np.arctan((ellipse.semi_x_axis / ellipse.semi_y_axis) * np.tan(_theta))


    # We only want quarter of the circumference...
    arc_length_q = ellipse.circumference / 4

    if from_t is None:
        arc_length_start = arc_length_q
        arc_length_end = 0.0
    else:
        arc_length_start = arc_length_q - calculate_ellipse_arc(ellipse, from_t)
        arc_length_end = 0.0

    equiarc_array_of_t = []

    # elliptic integrals work from the y-axis so start from the whole thing to get theta=0
    # for arc_length in np.linspace( arc_length_q, 0.0, num ):
    for arc_length in np.linspace(arc_length_start, arc_length_end, num):
        # get the elliptic angle (from the y-axis to the point)
        phi = EllipticIntegralHelper.find_phi_for_arc_length(ellipse, arc_length)

        t = EllipticIntegralHelper.convert_phi_to_t(phi)

        equiarc_array_of_t.append(t)

    return equiarc_array_of_t


def calculate_pts_for_equi_arc_spaced_t(ellipse, num):
    """
    Calculate points on an ellipse quarter so that it is divided into arcs of equal length

    :param ellipse:
    :param num:
    :return:
    """

    # calculate the equally spaced parametric angles (wrt arc-length)
    t_values = calculate_equi_arc_spaced_t(ellipse, num)

    pts = [ellipse.calculate_ellipse_pt_for_t(t) for t in t_values]

    return pts


def calculate_pts_for_equi_spaced_t(ellipse, num):
    """
    Calculate points on an ellipse quarter such that the values of t are equally spaced
    :param ellipse object
    :param num:
    :return:
    """
    t_values = np.linspace(0, pi / 2, num)

    pts = [ellipse.calculate_ellipse_pt_for_t(t) for t in t_values]

    return pts


def calculate_polyline_offset_from_ellipse(ellipse, offset):
    """
    Create points on the (quarter) ellipse and calculate the 'offset' position using the normal
    The points on the ellipse are calculated based on curvature so that regions of high
    curvature are resolved satisfactorily

    :param ellipse:
    :param offset: from ellipse (can be negative)
    :return:
    """

    def calculate_curvature_range():
        """
        Calculate the min/max curvature with a curvature-dependent number
        of values in each subrange
        """
        max_num_per_range = 50

        num = 20

        # these will be max and min curvatures
        k_a, k_b = ellipse.calculate_curvature_for_t(0.0), ellipse.calculate_curvature_for_t(pi / 2)

        # space out buckets by curvature
        k_buckets = np.exp(np.linspace(np.log(k_a), np.log(k_b), num))

        k_range = np.array([])
        for k1, k2 in zip(k_buckets, k_buckets[1:]):
            # use more points for high curvature ranges
            num_per_range = min(max(abs(k1 - k2), num), max_num_per_range)

            # don't add the last one to avoid repetition
            k_range = np.append(k_range, np.linspace(k1, k2, num_per_range)[:-1])

        # need to add the very last one
        k_range = np.append(k_range, np.array([k_b, ]))

        return k_range

    def opt_fn(x, k):
        """ optimisation function """
        return abs(ellipse.calculate_curvature_for_t(x) - k)

    def calculate_t_values():
        """
        :return:
        """

        if ellipse.semi_x_axis == ellipse.semi_y_axis:
            t_values = np.linspace(0.0, pi / 2, 200)
        else:
            t_values = []

            for kappa in calculate_curvature_range():
                min_result = opt.minimize_scalar(opt_fn, args=(kappa,), bounds=(0.0, pi / 2), method='bounded')

                if min_result.success:
                    t_values.append(min_result.x)

            # only needed when b>a but here for safety
            t_values.sort()

            # make sure 0 and pi/2 are at the start and end of the range
            t_values = np.append([0.0, ], [ele for ele in t_values if 0.0 < ele < pi / 2])
            t_values = np.append(t_values, [pi / 2, ])

        return t_values

    offset_pts = []

    for t in calculate_t_values():
        pt = ellipse.calculate_ellipse_pt_for_t(t)
        n_v = ellipse.calculate_ellipse_normal_for_t(t)

        offset_pt = pt + offset * n_v
        offset_pts.append(offset_pt)

    return offset_pts


def calculate_ellipse_line_poi(ellipse, pt1, pt2):
    """
    Assumes points are in the x-y plane

    :param ellipse:
    :param pt1:
    :param pt2:
    :return: Point object
    :rtype: p.Point
    """

    maj_ax, min_ax = ellipse.semi_x_axis, ellipse.semi_y_axis

    # is pt1 on the ellipse?
    if abs((pt1.x / maj_ax) ** 2 + (pt1.y / min_ax) ** 2 - 1) < 1e-5:
        return pt1.x, pt1.y

    # is pt2 on the ellipse?
    if abs((pt2.x / maj_ax) ** 2 + (pt2.y / min_ax) ** 2 - 1) < 1e-5:
        return pt2.x, pt2.y

    a_sq = maj_ax ** 2
    b_sq = min_ax ** 2

    dx = pt2.x - pt1.x
    dy = pt2.y - pt1.y

    # coefficients of quadratic equation
    c_a = a_sq * dy * dy + b_sq * dx * dx
    c_b = 2 * (dx * pt1.x * b_sq + dy * pt1.y * a_sq)
    c_c = b_sq * pt1.x ** 2 + a_sq * pt1.y ** 2 - a_sq * b_sq

    l1, l2 = geom_utils.solve_quadratic(c_a, c_b, c_c)

    if (l1 < 0.0 and l2 < 0.0) or (l1 > 1.0 and l2 > 1.0):
        raise ValueError('The line joining {} and {} does not intersect with ellipse (maj={}, min={})'
                         .format(pt1, pt2, maj_ax, min_ax))

    def dist_from_range(_x, _a, _b):
        """ minimum difference between x and a, and x and b """
        return min(abs(_x - _a), abs(_x - _b))

    if 0 < l1 < 1:
        _l = l1
    elif 0 < l2 < 1:
        _l = l2
    else:
        # take the one nearest to the range
        _l = l1 if dist_from_range(l1, 0.0, 1.0) < dist_from_range(l2, 0.0, 1.0) else l2

    return pt1 + _l * (pt2 - pt1)


# Plotting

def plot_initialise(add_grid=True):
    """
    Initialise the plot

    :param add_grid:
    :return:
    """

    import matplotlib.pyplot as plt

    ax = plt.figure().add_subplot(111, aspect='equal')
    ax.set_xlabel('x')
    ax.set_ylabel('y')

    plt.grid(add_grid)

    return plt, ax


def plot_add_pts(plt, pts, point_type, label=None):
    """
    Add points to a plot

    :param plt:
    :param pts:
    :param point_type:
    :param label:
    :return:
    """

    h_plt = plt.plot([pt.x for pt in pts], [pt.y for pt in pts], point_type, label=label)

    return h_plt


def plot_show(plt, ax):
    """
    Show a plot
    :param plt:
    :param ax:
    :return:
    """

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels)

    plt.show()

    return


def plot_arc_vs_polar_angle():
    """
    Comparison plot
    :return:
    """

    major, minor, nn = 2.0, 1.0, 10

    ellipse = Ellipse(major, minor)

    plt, ax = plot_initialise()

    theta_pts = calculate_pts_for_equi_spaced_t(ellipse, nn)
    plot_add_pts(plt, theta_pts, '.k', 'equi-t')

    arc_pts = calculate_pts_for_equi_arc_spaced_t(ellipse, nn)
    plot_add_pts(plt, arc_pts, '+r', 'equi-arc')

    plt.xlim((-0.1, major * 1.1))
    plt.ylim((-0.1, minor * 1.1))

    plot_show(plt, ax)


def plot_ellipse_and_offset():
    """
    Plot ellipse
    :return:
    """

    aa, bb, nn, offset = 2.0, 0.5, 100, 1.0

    ellipse = Ellipse(aa, bb)

    plt, ax = plot_initialise()

    # plt.xlim( (-a / 20, a * 1.05 ) )
    # plt.ylim( (-b / 10, b * 1.1 ) )

    arc_pts = calculate_pts_for_equi_arc_spaced_t(ellipse, nn)
    plot_add_pts(plt, arc_pts, 'r', 'equi-arc')
    # plot_add_pts( plt, arc_pts, '+r' )

    offset_pts = calculate_polyline_offset_from_ellipse(ellipse, offset)
    plot_add_pts(plt, offset_pts, 'b', 'offset')
    # plot_add_pts( plt, offset_pts, '+b', 'offset' )

    plot_show(plt, ax)

    return


def test_ellipse_3d():
    """
    Test the Ellipse3D code
    :return:
    """

    # TODO: move to unit tests

    pt_c = p.Point(1, 1, 1)
    pt_1 = p.Point(2, 0, 0) + pt_c
    pt_2 = p.Point(0, 0, 1) + pt_c

    e3d = Ellipse3D(pt_1, pt_2, pt_c)

    pt = e3d.calculate_ellipse_pt_for_t(0)
    assert pt == pt_1

    pt = e3d.calculate_ellipse_pt_for_t(pi / 2)
    assert pt == pt_2

    print(e3d)


def main():
    """
    :return:
    """

    if 1 == 1:
        test_ellipse_3d()
    else:
        plot_arc_vs_polar_angle()
        plot_ellipse_and_offset()

    return


if __name__ == '__main__':
    main()
