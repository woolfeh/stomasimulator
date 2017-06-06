from __future__ import print_function

import math

import numpy as np
from math import sqrt

import multiprocessing as mp

from . import point as p
from . import line as l


def solve_quadratic(aa, bb, cc):
    """
    Solve the quadratic equation: aa * x**2 + bb * x + cc = 0

    :param aa:
    :param bb:
    :param cc:
    :return: the two roots (-ve sqrt first)
    """

    # bring 'bb' out of the sqrt to help avoid an overflow situation
    sqrt_b2m4ac = bb * sqrt(1.0 - 4.0 * (aa / bb) * (cc / bb))

    return (-bb - sqrt_b2m4ac) / aa / 2, \
           (-bb + sqrt_b2m4ac) / aa / 2


def scalar_product(vec1, vec2):
    return sum(a * b for a, b in zip(vec1.xyz, vec2.xyz))


def scalar_triple_product(vec1, vec2, vec3):
    return scalar_product(vec1, vector_product(vec2, vec3))


def vector_product(vec1, vec2):
    x = vec1.y * vec2.z - vec1.z * vec2.y
    y = vec1.z * vec2.x - vec1.x * vec2.z
    z = vec1.x * vec2.y - vec1.y * vec2.x
    return p.Point(x, y, z)


def find_point_of_intersection(line, polyline):
    """
    Find point of intersection of line with polyline
    Straightforward O(n) algorithm because we have a line and a polyline
    :param line: a Line object
    :type line: l.Line
    :param polyline: list of Point objects
    :return: a Point object (most likely not in union(pl1,pl2)) or None
    """

    def plot_lines():
        import matplotlib.pyplot as plt

        plt.plot([p.x for p in line], [p.y for p in line], 'b')
        plt.plot([p.x for p in line], [p.y for p in line], 'ob')

        plt.plot([p.x for p in polyline], [p.y for p in polyline], 'r')
        plt.plot([p.x for p in polyline], [p.y for p in polyline], 'or')

        plt.show()

    poi = None
    if polyline is not None and len(polyline) > 0:
        # handy for debugging
        # plot_lines()

        # consider each consecutive pair of points
        for segment_pts in zip(polyline[:-1], polyline[1:]):
            segment = l.Line(segment_pts[0], segment_pts[1])

            # check the line segment if its extent overlaps
            # this will sometimes reduce the number that get checked
            if line.extent.overlaps(segment.extent):
                poi = line.find_point_of_intersection(segment)
                if poi is not None:
                    break

    return poi


def calc_unit_normal(p1, p2, p3, *pts):
    """ Calculate the unit normal from the given points - multiple points may be given
        (e.g. using the * prefix) but only the first three will be used.

        The normal vector's direction is defined by the vector product of v1=(p2-p1) and v2=(p3-p1).
        This means that the vector will point towards you when the points in the triangle p1-p2-p3
        are arranged in an anticlockwise sense.

        E.g. Let p1=(0,0,0), p2=(1,0,0), p3=(0,1,0)

               y

               |
               3                So normal will point in +ve z direction, i.e. (0,0,1)
               |
               |
               1____2__  x

        Mathematically, v1=p2-p1=(1,0,0), v2=p3-p1=(0,1,0) and v1 x v2 = (0,0,1).
    """

    # Need at least three so take first three
    assert isinstance(p1, p.Point) and isinstance(p2, p.Point) and isinstance(p3, p.Point)

    return vector_product(p2 - p1, p3 - p1).unit()


def rotate_about_x(pts, angle_degrees):
    """ Angle sense is as if looking towards +inf
        i.e. 90 degrees rotation takes (0,1,0) to (0,0,-1)

        NOTE: The id of each point is not preserved (it's treated as a new point) """

    sin_t = math.sin(math.pi / 180. * angle_degrees)
    cos_t = math.cos(math.pi / 180. * angle_degrees)

    rot_pts = (p.Point(1, 0, 0), p.Point(0, cos_t, sin_t), p.Point(0, -sin_t, cos_t))

    pp = []

    for pt_i in pts:
        new_triad = [scalar_product(rot_pt, pt_i) for rot_pt in rot_pts]
        new_p = p.from_triad(new_triad)
        pp.append(new_p)

    return pp


def calc_tetrahedron_volume(p1, p2, p3, p4):
    """
    Calculate the volume of a tetrahedron given its four vertices, each of which needs
    to be an instance of Point. Use the * prefix to pass a list of 4 points.

    :rtype : float
    """

    # use scalar triple product to get parallelopiped volume (tet. volume is 1/6 of that)
    stp = scalar_triple_product(p1 - p4, p2 - p4, p3 - p4)

    return abs(stp) / 6.0


def calculate_volume(pts_dict, facets):
    """
    Calculate the volume enclosed by the facets whose nodal coordinates are stored in pts_dict.
    Each facet is defined by a set of node ids and these are the keys to the map

    :param pts_dict Dictionary mapping point id to the Point object
    :param facets List of facets where each entry is a list of the node ids that comprise the facet
    :rtype : float The volume (always positive)
    """

    v, a = calculate_volume_and_area(pts_dict, facets)

    return v


def _calculate_facet_contrib(pts):
    n = calc_unit_normal(*pts)

    (n_gamma, idx) = n.max_coord()
    ap = AxesPermutations012[idx]

    pi_1 = calc_pi_1(pts, ap)

    d_area = pi_1 / n_gamma

    if idx == 0:
        # Max. normal component in the x-direction

        pi_alpha = calc_pi_alpha(pts, ap)
        pi_beta = calc_pi_beta(pts, ap)

        n_alpha = ap.alpha(n)
        n_beta = ap.beta(n)

        w = -scalar_product(n, pts[0])

        int_gamma_da = -sum(a * b for a, b in zip((n_alpha, n_beta, w), (pi_alpha, pi_beta, pi_1))) / (
            n_gamma * n_gamma)

        int_x_da = int_gamma_da

    elif idx == 1:
        # Max. normal component in the y-direction

        pi_beta = calc_pi_beta(pts, ap)
        int_beta_da = pi_beta / n_gamma
        int_x_da = int_beta_da

    else:
        # Max. normal component in the z-direction

        pi_alpha = calc_pi_alpha(pts, ap)
        int_alpha_da = pi_alpha / n_gamma
        int_x_da = int_alpha_da

    d_vol = n.x * int_x_da

    return d_vol, d_area


def calculate_volume_and_area(pts_dict, facets, num_procs=4):
    """
    Calculate the volume and surface area enclosed by the facets whose nodal coordinates are stored in pts_dict.
    Each facet is defined by a set of node ids and these are the keys to the map

    :param pts_dict: Dictionary mapping point id to the Point object
    :param facets: List of facets where each entry is a list of the node ids that comprise the facet
    :param num_procs:
    :rtype : 2-tuple The volume (always positive) and the area (signed)
    """

    try:
        # create a list where each item is a list of the points for each facet
        fpd = [[pts_dict[node_id] for node_id in f] for f in facets]
    except TypeError:
        print('Expected the points to be a dictionary (got {}) amd the facets to be list-like (got {})'.format(
            type(pts_dict), type(facets)))
        raise

    import platform
    if platform.system() == 'Windows':
        # Windows doesn't like 'Pool'
        vol_area = [_calculate_facet_contrib(pts) for pts in fpd]
    else:
        pool = mp.Pool(num_procs)

        # use map_async - let it work out chunk sizes
        # when it was tested against map, imap and imap_unordered it was as fast/faster
        vol_area = []
        asy = pool.map_async(_calculate_facet_contrib, fpd, callback=vol_area.extend)
        asy.wait()

        # clean up
        pool.close()
        pool.join()

    vol = sum(r[0] for r in vol_area)
    area = sum(r[1] for r in vol_area)

    # the face nodes may be defined s.t. the normal is pointing into the domain so return the absolute value
    #
    return abs(vol), area


def calculate_surface_area(pts_dict, facets):
    area = 0.0

    for facet in facets:
        facet_pts = [pts_dict.get(node_id) for node_id in facet]

        d_area = calculate_polygon_area(facet_pts)
        area += d_area

    return area


def calculate_polygon_area(pts_list):
    """
    Calculate the area enclosed by a list of points

    :param pts_list: list of Point objects (must have 3 or more entries) that are assumed to be coplanar
    :type pts_list: list of Point
    :return: returns the enclosed area (absolute value is returned)
    :rtype: float
    """

    if pts_list is None or len(pts_list) < 3:
        raise ValueError('The list of points must have at least 3 entries')

    def calculate_triangle_area(p1, p2, p3):
        """ Calculate the area of a triangle given its vertices, each of which needs to
            be an instance of Point. Use the * prefix to pass a list of 3 points. """

        return vector_product(p3 - p1, p3 - p2).l2_norm() / 2.0

    def calculate_polygon_area_impl():
        n = calc_unit_normal(*pts_list)
        (n_gamma, idx) = n.max_coord()

        ap = AxesPermutations012[idx]

        poly_area = calc_pi_1(pts_list, ap) / n_gamma

        return poly_area

    if len(pts_list) == 3:
        area = calculate_triangle_area(*pts_list)
    else:
        area = calculate_polygon_area_impl()

    return abs(area)


def calc_pi_1(pts, ap):
    lp = len(pts)

    assert (lp > 2)

    if pts[-1] == pts[0]:
        # closed loop so reduce the number of points (the loop treats 'lp' as the number of distinct vertices)
        lp -= 1

    pi_1 = 0.
    for i in xrange(lp):
        pt1 = pts[i]
        pt2 = pts[(i + 1) % lp]

        alpha1 = ap.alpha(pt1)
        alpha2 = ap.alpha(pt2)
        beta1 = ap.beta(pt1)
        beta2 = ap.beta(pt2)

        pi_1 += (beta2 - beta1) * (alpha2 + alpha1)

    return pi_1 / 2.


def calc_pi_alpha(pts, ap):
    lp = len(pts)

    assert (lp > 2)

    if pts[-1] == pts[0]:
        # closed loop so reduce the number of points (the loop treats 'lp' as the number of distinct vertices)
        lp -= 1

    pi_a = 0.
    for i in xrange(0, lp):
        pt1 = pts[i]
        pt2 = pts[(i + 1) % lp]

        alpha1 = ap.alpha(pt1)
        alpha2 = ap.alpha(pt2)
        beta1 = ap.beta(pt1)
        beta2 = ap.beta(pt2)

        pi_a += (beta2 - beta1) * (alpha2 * (alpha2 + alpha1) + alpha1 * alpha1)

    return pi_a / 6.


def calc_pi_beta(pts, ap):
    lp = len(pts)

    assert (lp > 2)

    if pts[-1] == pts[0]:
        # closed loop so reduce the number of points (the loop treats 'lp' as the number of distinct vertices)
        lp -= 1

    pi_b = 0.
    for i in xrange(0, lp):
        pt1 = pts[i]
        pt2 = pts[(i + 1) % lp]

        alpha1 = ap.alpha(pt1)
        alpha2 = ap.alpha(pt2)
        beta1 = ap.beta(pt1)
        beta2 = ap.beta(pt2)

        pi_b += (alpha2 - alpha1) * (beta2 * (beta2 + beta1) + beta1 * beta1)

    return -pi_b / 6.


class AxesPermute(object):
    """
    Permute the Cartesian axes to the alpha-beta-gamma parametric ones suitable for volume and
    surface area calculations.
    """

    def __init__(self, xyz_axis_char='x'):
        #                     x  y  z    y  z  x    z  x  y
        self.index_order = ((0, 1, 2), (1, 2, 0), (2, 0, 1))

        self.index = 0
        self.set_alpha_axis(xyz_axis_char)

    def set_alpha_axis(self, xyz_axis_char):
        assert (xyz_axis_char == 'x' or xyz_axis_char == 'y' or xyz_axis_char == 'z')

        # Examples:
        #   (x,y,z) -> (x,y,z) = (a,b,g)
        #   (x,y,z) -> (y,z,x) = (a,b,g)
        #   (x,y,z) -> (z,x,y) = (a,b,g)

        self.index = ord(xyz_axis_char) - ord('x')

    def alpha(self, p):
        return p.xyz[self.index_order[self.index][0]]

    def beta(self, p):
        return p.xyz[self.index_order[self.index][1]]


AxesPermutations012 = {
    0: AxesPermute('y'),
    # Max. component of the normal is in the x-direction so permute axes (x,y,z) -> (y,z,x) = (a,b,g)
    1: AxesPermute('z'),
    # Max. component of the normal is in the y-direction so permute axes (x,y,z) -> (z,x,y) = (a,b,g)
    2: AxesPermute('x')
    # Max. component of the normal is in the z-direction so permute axes (x,y,z) -> (x,y,z) = (a,b,g)
}


def calculate_euler_characteristic(surface):
    """
    Calculate the Euler characteristic for a polygonal surface, i.e. chi = V - E + F
    where
      V is the number of vertices
      E is the number of edges, and
      F is the number of faces

    The input, surface, is a list of faces that make up the surface.
    For closed surfaces chi == 2

    :rtype tuple of (chi, node set, edge set)
    """

    class MinMaxPair(object):
        def __init__(self, a, b):
            self.a = min(a, b)
            self.b = max(a, b)

        def __eq__(self, other):
            return isinstance(other, type(self)) and (self.a == other.a) and (self.b == other.b)

        def __hash__(self):
            return hash((self.a, self.b))

        def __repr__(self):
            return "({}, {})".format(self.a, self.b)

        @property
        def id1(self):
            return self.a

        @property
        def id2(self):
            return self.b

    node_set = set()
    edge_set = set()

    for facet in surface:
        for node_id in facet:
            node_set.add(node_id)

        lf = len(facet)
        for i in xrange(lf):
            edge_set.add(MinMaxPair(facet[i], facet[(i + 1) % lf]))

    num_nodes = len(node_set)
    num_edges = len(edge_set)
    num_faces = len(surface)

    chi = num_nodes - num_edges + num_faces

    return chi, node_set, edge_set


def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians (uses the Euler-Rodrigues formula).

    :param axis: 3-tuple
    :param theta: angle in radians
    :return: 3x3 rotation matrix
    """

    axis = np.asarray(axis)
    theta = np.asarray(theta)

    axis = axis / sqrt(np.dot(axis, axis))

    a = math.cos(theta / 2.0)

    b, c, d = -axis * math.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d

    rot_mat = np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                        [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                        [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])

    for item in np.nditer(rot_mat, op_flags=['readwrite']):
        item *= 0.0 if abs(item) < 1e-10 else 1.0

    return rot_mat


def rotate_3d(v, rot_matrix):
    rot_v = np.dot(rot_matrix, v)
    return rot_v


def rotation_interactive():
    import numpy as np

    v = [3, 5, 0]
    # axis = [4, 4, 1]
    axis = [0, 0, 1]
    theta = np.pi / 2

    r = rotation_matrix(axis, theta)
    v_r = rotate_3d(v, r)

    print('Rotation matrix is {}'.format(r))
    print('Vector is {}'.format(v))
    print('Rotated vector is {}'.format(v_r))


def lines_interactive():
    l1 = l.Line(p.Point(0, 0, 0), p.Point(1, 0, 0))

    # l2 = l.Line( p.Point(0, -0.5, 0), p.Point(1, 0.5, 0) )
    # l2 = l.Line( p.Point(0, -0.25, 0), p.Point(1, -0.25, 0) )
    l2 = l.Line(p.Point(0, -0.25, 0), p.Point(1, 0.1, 0))

    poi1 = l1.find_point_of_intersection(l2)
    print(poi1)

    l3 = l.Line(p.Point(0, 1, 0), p.Point(1, 1, 0))

    poi2 = l1.find_point_of_intersection(l3)
    print(poi2)

    polyline = (l1.pt1, l1.pt2, l3.pt1, l3.pt2)

    poi3 = find_point_of_intersection(l2, polyline)
    print(poi3)

    import matplotlib.pyplot as plt
    for pts in zip(polyline[:-1], polyline[1:]):
        plt.plot([pt.x for pt in pts], [pt.y for pt in pts], '-r')

    plt.plot([pt.x for pt in l2], [pt.y for pt in l2], '-g')

    for pt in (poi1, poi2, poi3):
        if pt is not None:
            plt.plot(pt.x, pt.y, 'o')

    plt.xlim([-0.5, 1.5])
    plt.ylim([-0.5, 1.5])

    plt.show()


def main():
    rotation_interactive()
    lines_interactive()


if __name__ == '__main__':
    main()
