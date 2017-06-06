####################################################################################################

from stomasimulator.geom.point import Point

import stomasimulator.mesh.elements as elem

####################################################################################################

class TestElements(object):

    def get_hex8(self, a=1.0, b=1.0, c=1.0, skew=0.0):

        pts = ( Point(0, 0, skew, allocate_id=True ),
                Point(a, 0, 0, allocate_id=True ),
                Point(a, b, 0, allocate_id=True ),
                Point(0, b, 0, allocate_id=True ),
                Point(0, 0, c + skew, allocate_id=True ),
                Point(a, 0, c + skew, allocate_id=True ),
                Point(a, b, c + skew, allocate_id=True ),
                Point(0, b, c + skew, allocate_id=True ) )


        pts_map = { pt.id: pt for pt in pts }

        ele = elem.Hex8Element( [ pt.id for pt in pts ] )

        return ele, pts_map

    def test_jacobian_cube(self):
        ele, pts_map = self.get_hex8()

        Jratio = ele.jacobian_ratio( pts_map )

        assert Jratio == 1.0

    def test_jacobian_rect(self):
        ele, pts_map = self.get_hex8( a=1.0, b=2.0, c=3.0 )

        Jratio = ele.jacobian_ratio( pts_map )

        assert Jratio == 1.0

    def test_jacobian_rect_skew(self):
        ele, pts_map = self.get_hex8( skew=1.0 )

        # tweak one point to make the faces non-parallel
        pts_map[ 10 ] = Point( 0, 0, -1 )

        Jratio = ele.jacobian_ratio( pts_map )

        assert Jratio != 1.0

####################################################################################################
