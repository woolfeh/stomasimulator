
from . import point as p


class BoxExtent( object ):
    """
    'Extent' represents a box in 3D Cartesian space with sides parallel to the xy, xz and yz planes.
    The box can be constructed from two 'Point's or from a list of 'Point's.

    Instance variables: p1 and p2 (p1 has the min. x, y and z coordinate)

    The extent can be queried to see if a point lies within the extent.
    """

    def __init__( self, pt1, pt2 ):
        self.p1 = p.Point( min( pt1.x, pt2.x ), min( pt1.y, pt2.y ), min( pt1.z, pt2.z ) )
        self.p2 = p.Point( max( pt1.x, pt2.x ), max( pt1.y, pt2.y ), max( pt1.z, pt2.z ) )

    def __str__( self ):
        return 'Extent: x in [{:8.4f}, {:8.4f}], y in [{:8.4f}, {:8.4f}], z in [{:8.4f}, {:8.4f}], (dx,dy,dz)=( {:8.4f}, {:8.4f}, {:8.4f} )'\
            .format( self.min_x, self.max_x, self.min_y, self.max_y, self.min_z, self.max_z, self.dx, self.dy, self.dz )

    def __repr__( self ):
        return str( self )

    @property
    def min_x( self ):
        return self.p1.x

    @property
    def min_y( self ):
        return self.p1.y

    @property
    def min_z( self ):
        return self.p1.z

    @property
    def mid_x( self ):
        return 0.5 * ( self.min_x + self.max_x )

    @property
    def mid_y( self ):
        return 0.5 * ( self.min_y + self.max_y )

    @property
    def mid_z( self ):
        return 0.5 * ( self.min_z + self.max_z )

    @property
    def mid_pt( self ):
        return p.Point( self.mid_x, self.mid_y, self.mid_z )

    @property
    def max_x( self ):
        return self.p2.x

    @property
    def max_y( self ):
        return self.p2.y

    @property
    def max_z( self ):
        return self.p2.z

    @property
    def dx( self ):
        return abs( self.max_x - self.min_x )

    @property
    def dy( self ):
        return abs( self.max_y - self.min_y )

    @property
    def dz( self ):
        return abs( self.max_z - self.min_z )

    @property
    def mid_point( self ):
        return 0.5 * ( self.p1 + self.p2 )

    def contains_point( self, p ):
        return self.min_x <= p.x <= self.max_x and \
            self.min_y <= p.y <= self.max_y and \
            self.min_z <= p.z <= self.max_z

    def contains_extent(self, ext):
        return self.contains_point( ext.p1 ) and self.contains_point( ext.p2 )

    def overlaps(self, ext):
        return not ( ext.max_x < self.min_x or self.max_x < ext.min_x or
                     ext.max_y < self.min_y or self.max_y < ext.min_y or
                     ext.max_z < self.min_z or self.max_z < ext.min_z )

    def get_contained_points( self, pts ):
        subset = [ ]

        for pt in pts:
            if self.contains_point( pt ):
                subset.append( pt )

        return subset

    def grow( self, dx=0., dy=0., dz=0. ):
        dp = 0.5 * p.Point( abs( dx ), abs( dy ), abs( dz ) )
        self.p1 = self.p1 - dp
        self.p2 = self.p2 + dp

        return self

    @classmethod
    def from_points( cls, pts ):

        min_x = min( pts, key=lambda pt: pt.x ).x
        min_y = min( pts, key=lambda pt: pt.y ).y
        min_z = min( pts, key=lambda pt: pt.z ).z

        max_x = max( pts, key=lambda pt: pt.x ).x
        max_y = max( pts, key=lambda pt: pt.y ).y
        max_z = max( pts, key=lambda pt: pt.z ).z

        p1 = p.Point( min_x, min_y, min_z)
        p2 = p.Point( max_x, max_y, max_z)

        return cls( p1, p2 )

    @classmethod
    def from_point(cls, pt ):
        return cls( pt, pt )


class SphericalExtent( object ):
    """
    An 'Extent' that represents a sphere. Has a BoxExtent that contains the sphere to help.

    The extent is parameterised by the centre and a radius.

    When checking to see whether a point lies in (or on) the sphere the BoxExtent is checked first.
    Only if the point lies in the BoxExtent is the distance to the sphere's centre checked.
    """

    def __init__( self, centre, radius ):
        self.centre = centre
        self.radius = radius
        self.boxExtent = None

        self.reset_box_extent( )

    def reset_box_extent( self ):
        dp = p.Point( self.radius, self.radius, self.radius )
        p1 = self.centre - dp
        p2 = self.centre + dp
        self.boxExtent = BoxExtent( p1, p2 )

    def grow( self, dr ):
        if dr == 0.0:
            return self

        assert ( self.radius + dr > 0 )

        self.radius += dr
        self.reset_box_extent( )

        return self

    def get_contained_points( self, pts ):
        subset = [ ]

        for pt in pts:
            if self.contains_point( pt ):
                subset.append( pt )

        return subset

    def contains_point( self, pt ):
        if self.boxExtent.contains_point( pt ):
            d = self.centre.distance( pt )
            if d <= self.radius:
                return True

        return False
