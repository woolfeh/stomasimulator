""" Mesh configuration information """

from __future__ import print_function

from math import pi

class MeshConfig( object ):
    """ Mesh parameters for the stoma """

    def __init__( self, stoma_cfg ):
        """
        Set the number of points that will be used to discretise the stoma mesh.

        The number of slices for half a guard cell and the number of points on a quarter
        of the circumference are hard-coded. The number of points through the guard cell wall
        reflects the wall thickness.

        :param stoma_cfg:
        """

        self.stoma_cfg = stoma_cfg

        wall_thicknesses = ( stoma_cfg.wall_thickness.dorsal,
                             stoma_cfg.wall_thickness.ventral,
                             stoma_cfg.wall_thickness.periclinal )

        num_w = max( 4, int( round( 4 * max( wall_thicknesses ) ) ) )

        self.num_pts_on_qtr_circumference = 23
        self.num_pts_through_wall = num_w
        self.num_slices = 30

        self.num_pts_tip_wall_radius = self._num_tip_wall_radial_points()

    @property
    def num_pts_on_semi_circumference( self ):
        """ Number of points on half the circumference of a GC """
        return 2 * self.num_pts_on_qtr_circumference - 1

    @property
    def num_pts_on_epidermis_vert(self):
        """ Number of points on the vertical surface of the epidermal wall """

        # TODO use a parameter for the wall instead of num_pts_on_semi_circumference
        return self.num_pts_on_semi_circumference

    def _num_tip_wall_radial_points( self ):
        """
        Calculate the number of pts in the tip wall such that the elements are square(ish).

        half-circumference = pi * tip_radius = ds * num_elements_circ
        =>  ds = pi * tip_radius / num_elements_circ
        =>  num elements on tip radius = tip radius / ds

        :return:
        """

        scfg = self.stoma_cfg

        tip_length_inner = scfg.tip_length - scfg.wall_thickness.dorsal - scfg.wall_thickness.ventral
        tip_height_inner = scfg.tip_height - 2 * scfg.wall_thickness.periclinal

        min_tip_pm_radius = 0.5 * min( tip_length_inner, tip_height_inner )

        ds = pi * min_tip_pm_radius / ( self.num_pts_on_semi_circumference - 1)
        num_pts = max( 2, int( min_tip_pm_radius / ds ) + 1 )

        return num_pts

    def num_pts_through_tip_wall( self, inner_ventral_pts ):
        """
        :param inner_ventral_pts: list of Point objects that lie on the inside of the ventral wall
        :return:
        """

        th_ventral = self.stoma_cfg.wall_thickness.ventral
        th_tip = self.stoma_cfg.wall_thickness.polar

        # permit 10% margin in addition to the stated tip wall thickness
        thickness = 1.1 * th_tip

        # distance from origin to inside mid gc section (in equatorial plane)
        b_pm_v = self.stoma_cfg.pore_width / 2 + th_ventral

        if b_pm_v < thickness:
            # the tip wall is thicker than the PM aperture!
            num_pts = max( 2, self.num_slices / 4 )
        else:
            # count the number of inner ventral points that lie *below* the tip wall thickness
            pts_below_wall_thickness = ( vi_pt for vi_pt in inner_ventral_pts if vi_pt.y < th_tip )

            # a trick to perform len(.) on a generator
            num_pts = 1 + sum( 1 for _ in pts_below_wall_thickness )

        return num_pts

    def print_detail( self ):
        """ Output information about the meshing parameters """

        print( ' Mesh parameters' )
        print( ' ---------------' )
        print( '# slices:         per half GC : {}'.format( self.num_slices ) )
        print( '# pts:      through cell wall : {}'.format( self.num_pts_through_wall ) )
        print( '# pts: per half-circumference : {}'.format( self.num_pts_on_semi_circumference ) )
        print( ' ' )


if __name__ == '__main__':
    pass
