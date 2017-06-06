####################################################################################################

import math

import pytest

import stomasimulator.geom.point as p
import stomasimulator.geom.geom_utils as g

####################################################################################################

class TestMethods( object ):
    ## Test normal vector computations

    def test_scalar_product( self ):
        p1 = p.Point( 1, 2, 3 )
        p2 = p.Point( 2, 3, 4 )

        assert g.scalar_product( p1, p2 ) == 1. * 2. + 2. * 3. + 3. * 4.

    def test_vector_prod( self ):
        # Two vectors corresponding to edges of tetrahedron in the positive octant
        v1 = p.Point( 2, 0, -2 )
        v2 = p.Point( 0, 2, -2 )

        assert g.vector_product( v1, v2 ) == p.Point( 4, 4, 4 )

    def test_vector_prod_unit( self ):
        # Two vectors corresponding to edges of tetrahedron in the positive octant
        v1 = p.Point( 2, 0, -2 )
        v2 = p.Point( 0, 2, -2 )

        k = 4. * math.sqrt( 3. )
        assert g.vector_product( v1, v2 ).unit() == p.Point( 4. / k, 4. / k, 4. / k )

    # can't compute normal if points are collinear
    def test_collinear( self ):
        pts = [ p.Point( 0, 0, 0 ), p.Point( 1, 0, 0 ), p.Point( 2, 0, 0 ) ]

        with pytest.raises( ArithmeticError ):
            g.calc_unit_normal( *pts )

    # normal should be in x direction for points in y-z plane
    def test_normal_x( self ):
        pts = [ p.Point( 0, 2, 0 ), p.Point( 0, 0, 2 ), p.Point( 0, -1, -1 ) ]

        assert g.calc_unit_normal( *pts ) == p.Point( 1, 0, 0 )

    # normal should be in y direction for points in x-z plane
    def test_normal_y( self ):
        pts = [ p.Point( -2, 0, 0 ), p.Point( 0, 0, 4 ), p.Point( 1, 0, -1 ) ]

        assert g.calc_unit_normal( *pts ) == p.Point( 0, 1, 0 )

    # normal should be in z direction for points in x-y plane
    def test_normal_z( self ):
        pts = ( p.Point( 1, 1, 0 ), p.Point( 2, 1, 0 ), p.Point( 4, 4, 0 ) )

        assert g.calc_unit_normal( *pts ) == p.Point( 0, 0, 1 )

    ## Volumes and areas

    def test_tetrahedron_volume( self ):
        pts = ( p.Point( 1, 0, 0 ), p.Point( 0, 1, 0 ), p.Point( 0, 0, 0 ), p.Point( 0, 0, 1 ) )

        assert g.calc_tetrahedron_volume( *pts ) == 1. / 6

    def test_triangle_area( self ):
        pts = ( p.Point( 1, 0, 0 ), p.Point( 0, 1, 0 ), p.Point( 1, 1, 0 ) )

        assert g.calculate_polygon_area( pts ) == 0.5

    def test_triangle_area_1( self ):
        pts = ( p.Point( 1, 0, 0 ), p.Point( 0, 1, 0 ), p.Point( 0, 0, 1 ) )

        assert g.calculate_polygon_area( pts ) == math.sqrt( 3. ) / 2.

    def test_area_sq_open( self ):
        c = ( p.Point( 0, 0, 0 ), p.Point( 2, 0, 0 ), p.Point( 2, 2, 0 ), p.Point( 0, 2, 0 ) )

        assert g.calculate_polygon_area( c ) == 4.

    def test_area_sq_closed( self ):
        c = [ p.Point( -2, -2, 0 ), p.Point( -1, -2, 0 ), p.Point( -1, -1, 0 ), p.Point( -2, -1, 0 ) ]

        c.append( c[ 0 ] )

        assert g.calculate_polygon_area( c ) == 1.

    def test_area_tri_open( self ):
        c = [ p.Point( 1, 1, 0 ), p.Point( 2, 1, 0 ), p.Point( 2, 2, 0 ) ]

        assert g.calculate_polygon_area( c ) == 0.5

    def test_area_tri_closed( self ):
        c = [ p.Point( 1, 1, 0 ), p.Point( 2, 1, 0 ), p.Point( 2, 2, 0 ) ]
        c.append( c[ 0 ] )

        assert g.calculate_polygon_area( c ) == 0.5

    def test_area_normal_towards_x( self ):
        c = [ p.Point( 3, -1, 0 ), p.Point( 3, 1, 0 ), p.Point( 0, 1, 4 ),
              p.Point( 0, -1, 4 ) ]

        assert g.calculate_polygon_area( c ) == 10.

    def test_area_normal_towards_y( self ):
        c = [ p.Point( 1, 3, 0 ), p.Point( -1, 3, 0 ), p.Point( -1, 0, 4 ),
              p.Point( 1, 0, 4 ) ]

        assert g.calculate_polygon_area( c ) == 10.

    def test_area_normal_towards_z( self ):
        c = [ p.Point( 3, -1, 0 ), p.Point( 3, 1, 0 ), p.Point( 0, 1, 1 ),
              p.Point( 0, -1, 1 ) ]

        assert g.calculate_polygon_area( c ) == 2. * math.sqrt( 10. )

    def test_calculate_euler_characteristic_1( self ):
        """
        Square     1--2
                   |  |
                   4--3
        """
        surface = ( (1, 2, 3, 4), )
        chi, d1, d2 = g.calculate_euler_characteristic( surface )

        assert chi == 1

    def test_calculate_euler_characteristic_2( self ):
        """
        Square     1--2--3
                   |  |  |
                   4--5--6
                   |  |  |
                   7--8--9
        """
        surface = ( (1, 2, 5, 4), (2, 3, 6, 5), (4, 5, 8, 7), (5, 6, 9, 8) )
        chi, d1, d2 = g.calculate_euler_characteristic( surface )

        assert chi == 1

    def test_calculate_euler_characteristic_3( self ):
        """
        Cube (hidden vertex is 5)
                  8--7
                 /  /|
                4--3 6
                |  |/
                1--2
        """
        surface = (
        (1, 2, 3, 4), (5, 6, 7, 8), (2, 6, 7, 3), (1, 5, 8, 4), (1, 2, 6, 5), (4, 3, 7, 8) )
        chi, d1, d2 = g.calculate_euler_characteristic( surface )

        assert chi == 2

    def test_calculate_euler_characteristic_4( self ):
        """
        Hole-y cube, i.e. topologically equivalent to a cylinder
        Hidden vertex is 5 - (1,2,3,4) and (5,6,7,8) are not faces)
                  8--7
                 /  /|
                4--3 6
                |  |/
                1--2
        """
        surface = ( (2, 6, 7, 3), (1, 5, 8, 4), (1, 2, 6, 5), (4, 3, 7, 8) )
        chi, d1, d2 = g.calculate_euler_characteristic( surface )

        assert chi == 0

    # v. simple volume calculations

    def test_vol_tetrahedron( self ):
        p1 = p.Point( 1, 0, 0 )
        p2 = p.Point( 0, 1, 0 )
        p3 = p.Point( 0, 0, 0 )
        p4 = p.Point( 0, 0, 1 )

        volume = g.calc_tetrahedron_volume( p1, p2, p3, p4 )

        assert volume == 1.0 / 6.0

    def test_vol_tetrahedron_1( self ):
        # Volume of a tetrahedron but using the general faces and vertices formula

        p1 = p.Point( 1, 0, 0, 1 )
        p2 = p.Point( 0, 1, 0, 2 )
        p3 = p.Point( 0, 0, 0, 3 )
        p4 = p.Point( 0, 0, 1, 4 )

        vertices = { 1: p1, 2: p2, 3: p3, 4: p4 }

        # Need a list of faces: vertices are labelled 1,2,3,4 so faces are 1-2-4, 2-3-4, 1-4-3 and 1-3-2
        # (all traversed anti-clockwise when looking from outside).
        #
        faces = ( (1, 2, 4), (2, 3, 4), (1, 4, 3), (1, 3, 2) )

        volume = g.calculate_volume( vertices, faces )

        assert volume == 1.0 / 6.0


    def test_vol_irregular_cuboid_1( self ):
        # Calculate volume of irregular rectangular cuboid using projection method
        # - node indices are as per FEBio representation of 8-node element

        a, b, c = 2., 3., 4.

        vertices, faces = self.get_cuboid( a, b, c )

        vertices_dict = dict( )
        for vertex in vertices:
            vertices_dict[ vertex.id ] = vertex

        vol = g.calculate_volume( vertices_dict, faces )

        assert vol == a * b * c

    def test_vol_irregular_cuboid_2( self ):
        # Calculate volume of irregular rectangular cuboid using projection method
        # - node indices are as per FEBio representation of 8-node element

        a, b, c = 2., 3., 4.

        vertices, faces = self.get_cuboid( a, b, c )

        vv = g.rotate_about_x( vertices, 45 )

        # id is not preserved so use the one from above
        vertex_dict = dict( )
        for idx in range( len( vertices ) ):
            vertex_dict[ vertices[ idx ].id ] = vv[ idx ]

        vol = g.calculate_volume( vertex_dict, faces )

        assert vol == a * b * c

    def get_cuboid( self, a, b, c ):

        # node indices are as per FEBio representation of 8-node element

        vertices = list()
        vertices.append( p.Point( a, 0, 0 ) )
        vertices.append( p.Point( a, b, 0 ) )
        vertices.append( p.Point( 0, b, 0 ) )
        vertices.append( p.Point( 0, 0, 0 ) )
        vertices.append( p.Point( a, 0, c ) )
        vertices.append( p.Point( a, b, c ) )
        vertices.append( p.Point( 0, b, c ) )
        vertices.append( p.Point( 0, 0, c ) )

        # rely on the id so set it
        p.Point.point_id = 1

        for pt in vertices:
            pt.allocate_id()

        # List of faces
        faces = ( (1, 4, 3, 2), (5, 6, 7, 8), (1, 5, 8, 4), (2, 3, 7, 6), (3, 4, 8, 7), (1, 2, 6, 5) )

        return vertices, faces

####################################################################################################

    class TestRotationMatrix(object):
        def test_rotate_around_x(self):
            p = [ 1, 2, 3 ]
            rm = g.rotation_matrix( [1, 0, 0], math.pi / 2 )

            rot_p = g.rotate_3d( p, rm )

            assert rot_p[0] == p[0] and rot_p[1] == -p[2] and rot_p[2] == p[1]

        def test_rotate_around_y(self):
            p = [ 1, 2, 3 ]
            rm = g.rotation_matrix( [0, 1, 0], math.pi / 2 )

            rot_p = g.rotate_3d( p, rm )

            assert rot_p[0] == p[2] and rot_p[1] == p[1] and rot_p[2] == -p[0]

        def test_rotate_around_z(self):
            p = [ 1, 2, 3 ]
            rm = g.rotation_matrix( [0, 0, 1], math.pi / 2 )

            rot_p = g.rotate_3d( p, rm )

            assert rot_p[0] == -p[1] and rot_p[1] == p[0] and rot_p[2] == p[2]

        def test_rotate_around_x_to_xy_plane(self):
            import numpy as np

            rr = np.linspace( -4, 4, 17 )

            for x in rr:
                for y in rr:
                    for z in rr:
                        xyz = [ x, y, z ]

                        theta = np.arctan2( z, y )
                        rm = g.rotation_matrix( [1, 0, 0], -theta )

                        rot_xyz = g.rotate_3d( xyz, rm )

                        assert ( rot_xyz[1] >= 0.0 or abs( rot_xyz[1] ) < 1e-10 ) and abs( rot_xyz[2] ) < 1e-10

        def test_rotate_around_y_to_xy_plane(self):
            import numpy as np

            rr = np.linspace( -4, 4, 17 )

            for x in rr:
                for y in rr:
                    for z in rr:
                        xyz = [ x, y, z ]

                        theta = np.arctan2( x, z )
                        rm = g.rotation_matrix( [0, 1, 0], np.pi / 2 - theta )

                        rot_xyz = g.rotate_3d( xyz, rm )

                        assert ( rot_xyz[0] >= 0.0 or abs ( rot_xyz[0] ) < 1e-10 ) and abs( rot_xyz[2] ) < 1e-10

        def test_rotate_around_z_to_xz_plane(self):
            import numpy as np

            rr = np.linspace( -4, 4, 17 )

            for x in rr:
                for y in rr:
                    for z in rr:
                        xyz = [ x, y, z ]

                        theta = np.arctan2( y, x )
                        rm = g.rotation_matrix( [0, 0, 1], -theta )

                        rot_xyz = g.rotate_3d( xyz, rm )

                        assert ( rot_xyz[0] >= 0.0 or abs( rot_xyz[0] ) < 1e-10 ) and abs( rot_xyz[1] ) < 1e-10

####################################################################################################
