
import numpy as np
import numpy.linalg as la

from . import extent as ext


class Line( object ):
    """
    A 'Line' is defined by providing two points on the line. The line is stored
    parametrically by
      P(t) = A + t * ( B - A )
    where P are the points on the line, A and B are the two points defining the line, and
      t is the (scalar) parameter that takes A to B, with t in [0, 1] for points on the line
    """

    def __init__( self, pt1, pt2 ):
        self.pt1 = pt1
        self.pt2 = pt2
        self.extent = ext.BoxExtent( pt1, pt2 )

    def get_point( self, t ):
        return self.pt1 + t * ( self.pt2 - self.pt1 )

    def find_point_of_intersection(self, l):
        """
        Line 1: A + s(B-A)
        Line 2: C + t(D-C)

        [ B-A C-D ] [ s t ]' = [ C-A ]

        where B-A, C-D and C-A are column vectors

        s and t are the parameters for the respective lines.
        """
        assert isinstance( l , Line ), 'Parameter must be a Line object'

        def in_range( value, lb, ub, d=np.spacing( np.single(1) ) ):
            return lb - d <= value <= ub + d

        # Line 1 is A + t(B-A)
        A = np.array( self.pt1.xyz )
        B = np.array( self.pt2.xyz )

        # Line 2 is C + t(D-C)
        C = np.array( l.pt1.xyz )
        D = np.array( l.pt2.xyz )

        # B-A and C-D are row vectors here so transpose them...
        M = np.matrix( [ B-A, C-D ] ).transpose()
        RHS = C-A

        try:
            rtc = la.lstsq( M, RHS )
            X = rtc[0]
        except la.LinAlgError:
            # Hasn't converged so assume lines don't converge => return None
            return None

        the_pt = None

        # X[i] is line parameter for line i+1
        s = X[0]
        t = X[1]
        if in_range( s, 0.0, 1.0 ) and in_range( t, 0.0, 1.0 ):
            pt0 = self.get_point( s )
            pt1 = l.get_point( t )

            the_pt = pt0 if pt0.distance( pt1 ) < 1e-10 else None

        return the_pt

    def __iter__(self):
        yield self.pt1
        yield self.pt2
