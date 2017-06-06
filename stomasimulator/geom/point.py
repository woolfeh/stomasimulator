import scipy.spatial.distance as spd
import numpy as np

DEFAULT_ALLOCATE_ID = False


def from_triad(tri, allocate_id=None):
    assert isinstance(tri, (list, tuple, np.ndarray)) and len(tri) == 3
    return Point(tri[0], tri[1], tri[2], allocate_id=allocate_id)


class Point(object):
    """
    The class 'Point' represents a point in 3D Cartesian space.
    Standard operations are implemented, e.g. +, -, * (scalar product), norm
    Instance attributes: x y z
    """

    point_id = 1

    def __init__(self, x, y, z, allocate_id=False):
        """
        Create a 3D point

        :param x: x coordinate
        :param y: y coordinate
        :param z: z coordinate
        :param allocate_id: whether to assign an id to the point
        :return:
        """

        self._xyz = None
        self.xyz = (x, y, z)

        # TODO remove id from the class - it is context specific - Point should be just POD
        self.id = None
        if allocate_id is not None and allocate_id:
            self.allocate_id()

        return

    def is_id_set(self):
        return self.id is not None

    def allocate_id(self):
        if not self.is_id_set():
            self.id = Point.point_id
            Point.point_id += 1

        return self.id

    @property
    def xyz(self):
        return self._xyz

    @xyz.setter
    def xyz(self, new_xyz):
        assert len(new_xyz) == 3

        def clean_num(n):
            n = float(n)
            return 0.0 if abs(n) < 1e-10 else n

        self._xyz = [clean_num(_) for _ in new_xyz]

    @property
    def x(self):
        return self._xyz[0]

    @x.setter
    def x(self, xx):
        self._xyz[0] = xx

    @property
    def y(self):
        return self._xyz[1]

    @y.setter
    def y(self, yy):
        self._xyz[1] = yy

    @property
    def z(self):
        return self._xyz[2]

    @z.setter
    def z(self, zz):
        self._xyz[2] = zz

    # special methods

    def __eq__(self, other_pt):
        # 'id' is not checked - equality is purely spatial
        return isinstance(other_pt, Point) and np.all(np.equal(self.xyz, other_pt.xyz))

    __hash__ = None

    def __ne__(self, other):
        return not self == other

    def __str__(self):
        if self.id is None or self.id == 0:
            return '({}, {}, {})'.format(self.x, self.y, self.z)
        else:
            return '({}, {}, {})-id={}'.format(self.x, self.y, self.z, self.id)

    def __repr__(self):
        return self.__str__()

    def __add__(self, pt):
        return Point(self.x + pt.x, self.y + pt.y, self.z + pt.z)

    def __radd__(self, pt):
        if isinstance(pt, Point):
            return Point(self.x + pt.x, self.y + pt.y, self.z + pt.z)
        else:
            return Point(self.x + pt, self.y + pt, self.z + pt)

    def __sub__(self, pt):
        if isinstance(pt, Point):
            return Point(self.x - pt.x, self.y - pt.y, self.z - pt.z)
        raise TypeError("The operation '{} - {}' is undefined".format(self, pt))

    def __rsub__(self, pt):
        if isinstance(pt, Point):
            return Point(pt.x - self.x, pt.y - self.y, pt.z - self.z)
        raise TypeError("The operation '{} - {}' is undefined".format(pt, self))

    def __mul__(self, value):
        # multiply by a scalar
        if isinstance(value, (int, long, float)):
            return Point(value * self.x, value * self.y, value * self.z)

        return np.dot(self.xyz, value.xyz)

    def __rmul__(self, a):
        return self * a

    def __div__(self, a):
        return Point(self.x / a, self.y / a, self.z / a)

    def __neg__(self):
        return Point(-self.x, -self.y, -self.z)

    # user-defined methods

    def l2_norm(self):
        return spd.euclidean(self.xyz, (0, 0, 0))

    def distance(self, pt):
        return spd.euclidean(self.xyz, pt.xyz)

    def unit(self):
        mag = self.l2_norm()
        if mag == 0.0:
            raise ArithmeticError('Magnitude is zero - unable to calculate unit vector')

        return Point(self.x / mag, self.y / mag, self.z / mag)

    def max_coord(self):
        # Will return (value, index-base-0) for the maximum of the x, y and z values

        (mag, i) = max((abs(v), i) for i, v in enumerate(self.xyz))

        return self.xyz[i], i


x_axis = Point(1, 0, 0)
y_axis = Point(0, 1, 0)
z_axis = Point(0, 0, 1)
xyz_axes = (x_axis, y_axis, z_axis)

xyz_axes_dict = dict(x=x_axis, y=y_axis, z=z_axis)


def xyz_axis(i):
    assert i in (0, 1, 2)
    return xyz_axes[i]


def translate_pt_by_xyz(pt, xyz):
    """
    Move the point - creates a new point

    :param pt: the point to move
    :param xyz: by this 3-tuple
    :return: a new point
    """

    assert isinstance(pt, Point) and (isinstance(xyz, (list, tuple)) and len(xyz) == 3)

    xyz = [pt.xyz[i] + xyz[i] for i in range(3)]

    return Point(xyz[0], xyz[1], xyz[2])
