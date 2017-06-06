####################################################################################################

from math import pi

import numpy as np

import stomasimulator.geom.ellipse as el
from stomasimulator.geom.point import Point

####################################################################################################

TOLERANCE = 1e-10

####################################################################################################

class TestEllipse(object):
    """
    Test the code in ellipse.py

    Here, a and b are the semi- and semi-minor axes (x is assumed to be the major axis and so a>=b)
    """


    def test_calculate_m_for_ellipe(self):
        a, b = 3.0, 2.0

        ellipse = el.Ellipse( a, b )

        assert ellipse.m == 1.0 - ( b / a )**2

    def test_calculate_ellipse_x(self):
        a, b = 3.0, 2.0

        ellipse = el.Ellipse( a, b )

        y = b / 2
        x = ellipse.calculate_ellipse_x( y )

        assert abs( (x/a)**2 + (y/b)**2 - 1.0 ) < TOLERANCE

    def test_calculate_ellipse_y(self):
        a, b = 3.0, 2.0

        ellipse = el.Ellipse( a, b )

        x = a / 2
        y = ellipse.calculate_ellipse_y( x )

        assert abs( (x/a)**2 + (y/b)**2 - 1.0 ) < TOLERANCE

    def test_calculate_ellipse_xy_for_t(self):
        a, b = 3.0, 2.0

        ellipse = el.Ellipse( a, b )

        for t in ( 0.0, pi/4, pi/2 ):
            x, y = ellipse.calculate_ellipse_xy_for_t( t )
            assert abs( (x/a)**2 + (y/b)**2 - 1.0 ) < TOLERANCE

    def test_calculate_equi_arc_spaced_t(self):
        a, b = 2.0, 1.0

        ellipse = el.Ellipse( a, b )

        num = 100
        thetas = el.calculate_equi_arc_spaced_t( ellipse, num )

        # calculate the arc for each parametric angle
        cumulative_arcs = np.array( [ el.calculate_ellipse_arc( ellipse, theta ) for theta in thetas ] )

        # calculate the length of each arc
        arcs = cumulative_arcs[1:] - cumulative_arcs[:-1]

        # the arcs should all be the same size
        for idx in range(num-2):
            assert abs( arcs[idx] - arcs[idx+1] ) < TOLERANCE

    def test_find_phi_for_arc_length(self):
        # Test on a circle (the function is used by calculate_equi_arc_spaced_t so it's tested already)
        # This is just for completeness (and also makes sure circles work too!)

        # circle radius
        a = 2.0

        ellipse = el.Ellipse( a, a )

        # arc length at phi = 0, pi/4 and pi/2
        arc_lengths = [ 0.0, pi * a / 4, pi * a / 2 ]

        # get the elliptic angles, phi, that are measured from the y-axis
        phis = [ el.EllipticIntegralHelper.find_phi_for_arc_length( ellipse, l) for l in arc_lengths ]

        for i, arc_length in enumerate( arc_lengths ):
            assert arc_length == phis[i] * a

####################################################################################################

class TestEllipse3D( object ):

    def test_identity(self):
        a, b = 2, 1
        ellipse = el.Ellipse( a, b )
        c = Point( 0, 0, 0 )
        ellipse3D = el.Ellipse3D( Point( a, 0, 0), Point( 0, b, 0 ), c )

        assert ellipse.semi_x_axis == ellipse3D.axis_pt_1.x - c.x and \
               ellipse.semi_y_axis == ellipse3D.axis_pt_2.y - c.y

    def test_quadrant1_ellipse(self):
        a, b = 2, 1
        ellipse = el.Ellipse( a, b )

        c = Point( 1, 2, 3 )
        ellipse3D = el.Ellipse3D( Point( a, 0, 0) + c, Point( 0, b, 0 ) + c, c )

        assert ellipse.semi_x_axis == ellipse3D.axis_pt_1.x - c.x and \
               ellipse.semi_y_axis == ellipse3D.axis_pt_2.y - c.y

####################################################################################################
