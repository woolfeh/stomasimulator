####################################################################################################

import math
import pytest

import stomasimulator.geom.point as p

####################################################################################################

class TestPoint( object ):
    def test_pt_eq( self ):
        p1 = p.Point( 1, 2, 3 )
        p2 = p.Point( 1, 2, 3 )

        assert p1 == p2

    def test_pt_add( self ):
        p1 = p.Point( 1, 2, 3 )
        p2 = p.Point( 2, 3, 4 )

        assert p1 + p2 == p.Point( 3, 5, 7 )

    def test_pt_subtract( self ):
        p1 = p.Point( 1, 2, 3 )
        p2 = p.Point( 2, 4, 6 )

        assert p1 - p2 == p.Point( -1, -2, -3 )

    def test_pt_subtract_1( self ):
        pt = p.Point( 1, 2, 3 )

        with pytest.raises( TypeError ):
            pt - 1.

    def test_pt_subtract_2( self ):
        pt = p.Point( 1, 2, 3 )

        with pytest.raises( TypeError ):
            1. - pt

    def test_scalar_mult( self ):
        p1 = p.Point( 1, 2, 3 )
        k = 2.

        assert k * p1 == p.Point( k * p1.x, k * p1.y, k * p1.z )

    def test_scalar_mult2( self ):
        p1 = p.Point( 1, 2, 3 )
        k = 2.

        assert p1 * k == p.Point( k * p1.x, k * p1.y, k * p1.z )

    def test_div( self ):
        p1 = p.Point( 1, 2, 3 )
        k = 2.

        assert p1 / k == p.Point( 1 / k, 2 / k, 3 / k )

    def test_unary_minus( self ):
        p1 = p.Point( 1, 2, 3 )

        assert -p1 == p.Point( -1, -2, -3 )

    def test_l2_norm( self ):
        p1 = p.Point( 1, 2, 3 )

        assert p1.l2_norm( ) == math.sqrt( 1 + 4 + 9 )

    def test_distance( self ):
        p1 = p.Point( 1, 2, 3 )
        p2 = p.Point( 2, 3, 4 )

        assert p1.distance( p2 ) == math.sqrt( 3. )

    def test_xyz( self ):
        pt = p.Point( 1, 2, 3 )

        assert pt.xyz[0] == 1 and pt.xyz[1] == 2 and pt.xyz[2] == 3

    def test_max_coord_1( self ):
        pt = p.Point( 1, 2, 3 )

        assert pt.max_coord( ) == (3, 2)  # 2 because 0-based

    def test_max_coord_2( self ):
        pt = p.Point( 1, 3, 3 )

        assert pt.max_coord( ) == (3, 2)  # 2 because 0-based

    def test_max_coord_3( self ):
        pt = p.Point( -1, -4, 3 )

        assert pt.max_coord( ) == (-4, 1)  # 2 because 0-based

####################################################################################################
