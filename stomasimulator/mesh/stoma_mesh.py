""" Meshing specific to a stoma """

import numpy as np

from scipy import interpolate

from stomasimulator.geom import ( ellipse as el,
                                  geom_utils as g,
                                  line,
                                  point as p )

from stomasimulator.stomata import stoma_config as sc

from . import simple_mesh
from . import elements as elem

class StomaMesh(object):
    """ Representation of the mesh for a stoma - 'Mesh' object is an attribute (favour object composition) """

    def __init__(self, stoma_cfg):
        """
        Initialise stoma mesh object
        :param stoma_cfg: the stoma configuration object
        :type stoma_cfg: sc.StomaConfig
        """
        self.stoma_cfg = stoma_cfg

        # Generate the mesh...
        self.mesh = simple_mesh.Mesh( stoma_cfg )

    @property
    def mesh_cfg(self):
        """
        Short-cut to
        :return: the mesh configuration
        """
        return self.mesh.mesh_config

    def build(self, perform_mesh_checks):
        """
        Builds the mesh
        :param perform_mesh_checks:
        :return: the underlying mesh object
        :rtype: simple_mesh.Mesh
        """

        self._add_points( )
        self._add_elements( )

        if self.stoma_cfg.is_epidermal_wall_present:
            self._add_epidermal_cell_wall( )

        self._add_gc_reflections( )

        # set the boundary conditions - fix x at middle, y on polar walls and z for equatorial points
        self.mesh.set_bcs()

        # identify the facets on the inside and outside of the stoma - to apply the pressure etc.
        self.mesh.set_surface_facets( )

        if perform_mesh_checks:
            self.mesh.perform_simple_check( )

    def _add_points( self ):
        """ add the points for the guard cell quarter (in the octant with non-negative x, y and z) """

        self._add_equatorial_gc_points( )
        self._add_periclinal_gc_points( )
        self._add_intermediate_gc_points( )

        # the edges are added but are not used at present - the GmshWriter uses them
        # _add_gc_edges( mesh )

        if self.stoma_cfg.wall_thickness.polar > 0.0:
            self._add_tip_wall_points( )


    def _add_elements( self ):
        """ Add the elements """

        self._add_gc_elements( )

        if self.stoma_cfg.wall_thickness.polar > 0.0:
            self._add_tip_wall_elements( )

    def _add_epidermal_cell_wall( self ):
        """ Add the epidermal cell wall

        :param mesh:
        :return
        """

        # add an offset to see the wall away from the dorsal wall
        testing_offset = 0.01

        # TODO add a wall_height parameter
        wall_height_factor = 1.5
        wall_height = wall_height_factor * self.stoma_cfg.mid_gc_height

        num_pts_z = self.mesh_cfg.num_pts_on_epidermis_vert

        semi_maj_ax = self.stoma_cfg.dorsal_ellipse.semi_x_axis + 0.01 + testing_offset
        semi_min_ax = self.stoma_cfg.dorsal_ellipse.semi_y_axis + 0.01 + testing_offset

        th_e = self.stoma_cfg.wall_thickness.epidermal

        ellipse_i = el.Ellipse( semi_maj_ax, semi_min_ax )
        ellipse_o = el.Ellipse( semi_maj_ax + th_e, semi_min_ax + th_e )

        nn = 2 * self.mesh_cfg.num_slices

        pts_i = el.calculate_pts_for_equi_arc_spaced_t( ellipse_i, nn )
        pts_o = el.calculate_pts_for_equi_arc_spaced_t( ellipse_o, nn )

        for c_idx, z in enumerate( np.linspace( 0.0, wall_height, num_pts_z ) ):
            for slice_idx, pt_io in enumerate( zip( pts_i, pts_o ) ):
                for t_idx in (0, 1):
                    pt = p.Point( pt_io[ t_idx ].x, pt_io[ t_idx ].y, z )

                    self.mesh.epidermal_nodes_grid[ t_idx, c_idx, slice_idx ] = pt
                    self.mesh.add_node( pt )

        self.mesh.add_special_node_id( 'epidermis', self.mesh.epidermal_nodes_grid[ 0, 0, 0 ].id )

        # index order is thickness, circumference and then slice
        #
        indices = ((0, 0, 0), (1, 0, 0), (1, 0, 1), (0, 0, 1),
                   (0, 1, 0), (1, 1, 0), (1, 1, 1), (0, 1, 1))

        for circumf_idx in range( num_pts_z - 1 ):
            for slice_idx in range( nn - 1 ):
                nid_indices = [ (ii[ 0 ], ii[ 1 ] + circumf_idx, ii[ 2 ] + slice_idx) for ii in indices ]

                nids = [ self.mesh.epidermal_nodes_grid[ t, c, s ].id for t, c, s in nid_indices ]

                ele = elem.Hex8Element( nids )
                self.mesh.add_element( ele )

        return

    def _add_gc_reflections( self ):
        """ GC quarter has been built - build the rest using reflections... """

        if self.stoma_cfg.num_gcs is None:
            return

        # reflect in z=0 plane
        self.mesh.add_xy_reflection( )

        # reflect in x=0 plane
        self.mesh.add_yz_reflection( )

        if self.stoma_cfg.num_gcs == 2:
            self.mesh.add_xz_reflection( )

    def _add_equatorial_gc_points( self ):
        r"""
        Control points on equator are:
        ==============================

           +------------------------                y
           | wall                   \
           +-----------------------  \              ^
                                   \  \             |
                                    \  \            +-> x
           +----------------         \  \
           | wall           \         \  \
           +---------------  \         \  \
                           \  \         \  \
                           l--m----p----r--t
                           |  |    |    |  |
                           k--n----o----q--s

        Points 'l' and 't' will not be collinear

        :param mesh:
        :return:
        """

        def add_ventral_dorsal_pts():
            """ Add the ventral and dorsal points - they lie on ellipses ) """

            for ellipse, circumf_idx in ( ( self.stoma_cfg.ventral_ellipse, 0 ),
                                          ( self.stoma_cfg.dorsal_ellipse, -1 ) ):
                outer_pts = el.calculate_pts_for_equi_arc_spaced_t( ellipse, self.mesh_cfg.num_slices )
                self._add_points_along_axis( 0, circumf_idx, outer_pts )

        def add_thickness_pts():
            """ Add ventral-dorsal points that lie in the wall (working from outside to inside) """

            ventral_pts = self.mesh.nodes_grid[ 0, 0, : ]
            dorsal_pts = self.mesh.nodes_grid[ 0, -1, : ]

            for slice_idx, vd_pt_pair in enumerate( zip( ventral_pts, dorsal_pts ) ):
                v_pt, d_pt = vd_pt_pair[ 0 ], vd_pt_pair[ 1 ]

                # line joining the ventral to the dorsal point - sections are constructed so that these points
                # define a vertical plane (i.e. plane's normal lies in the x-y plane)
                l_vd = line.Line( v_pt, d_pt )

                for outer_pt, inner_cv, circumf_idx in ( ( v_pt, self.stoma_cfg.inner_ventral_curve, 0 ),
                                                         ( d_pt, self.stoma_cfg.inner_dorsal_curve, -1 ) ):

                    inner_pt = g.find_point_of_intersection( l_vd, inner_cv )

                    thickness_pts = [ outer_pt, ]

                    # calculate the internal wall points (if required) going from the outside to the inside
                    if self.mesh_cfg.num_pts_through_wall > 2:
                        out_in_line = line.Line( outer_pt, inner_pt )

                        for line_t in np.linspace( 0.0, 1.0, self.mesh_cfg.num_pts_through_wall )[ 1:-1 ]:
                            th_pt = out_in_line.get_point( line_t )
                            thickness_pts.append( th_pt )

                    thickness_pts.append( inner_pt )

                    self._add_points_along_thickness( circumf_idx, slice_idx, thickness_pts )

        def store_special_points():
            """ Store points on the inside and outside of the GC in the middle.
                These can be used to find the inner points/facets later """

            for key, indices in ( ( 'pressure_node_ids', (-1,  0, -1) ),
                                  ( 'dorsal_node_ids',   ( 0, -1, -1) ) ,
                                  ( 'ventral_node_ids',  ( 0,  0, -1) ) ):
                th_idx, circumf_idx, slice_idx = indices
                pt = self.mesh.nodes_grid[ th_idx, circumf_idx, slice_idx ]
                self.mesh.add_special_node_id( key, pt.id )

        add_ventral_dorsal_pts()
        add_thickness_pts()
        store_special_points()

    def _add_periclinal_gc_points( self ):
        """
        Add the periclinal points (outside -> inside) to the nodes grid

        :param mesh:
        :return:
        """

        def get_smoothed_z():
            """ Use the parametric angle for the pt to smooth the transition (using the sin function) from the tip to
                the middle - this results in a smooth outer surface around the mid-point """
            _t = self.stoma_cfg.mid_ellipse.calculate_ellipse_t_for_pt( mid_pt )
            _z = tip_para_z + (mid_para_z - tip_para_z) * np.sin( _t )

            return _z

        # z-coordinates at the tip and at the mid-point
        tip_para_z = 0.5 * self.stoma_cfg.tip_height
        mid_para_z = 0.5 * self.stoma_cfg.mid_gc_height

        # get the ventral and dorsal points from the grid
        ventral_wall_pts = self.mesh.nodes_grid[ 0, 0, : ]
        dorsal_wall_pts = self.mesh.nodes_grid[ 0, -1, : ]

        periclinal_idx = self.mesh_cfg.num_pts_on_qtr_circumference - 1

        for slice_idx, vd_pair in enumerate( zip( ventral_wall_pts, dorsal_wall_pts ) ):
            # periclinal point will lie on the mid ellipse but must also be coplanar with the slice
            mid_pt = el.calculate_ellipse_line_poi( self.stoma_cfg.mid_ellipse, vd_pair[ 0 ], vd_pair[ 1 ] )

            z = get_smoothed_z()

            # add the points through the periclinal wall
            z_values = np.linspace( z,
                                    z - self.stoma_cfg.wall_thickness.periclinal,
                                    self.mesh_cfg.num_pts_through_wall )

            for z_idx, z_coord in enumerate( z_values ):
                periclinal_pt = p.Point( mid_pt.x, mid_pt.y, z_coord )

                self.mesh.add_node( periclinal_pt )
                self.mesh.nodes_grid[ z_idx, periclinal_idx, slice_idx ] = periclinal_pt

        return

    def _add_intermediate_gc_points( self ):
        """
        Calculate the points on a GC cross-section so that they are coplanar:
        1) First we calculate the ventral and dorsal points from the equi-arc spaced angles).
        2) Next we calculate the points on the (ventral and dorsal) PM ellipse and the central
           ellipse by finding the intersection point of the line joining the ventral and dorsal
           points and the relevant ellipse

        By coplanar we mean in the same plane as that defined by line joining the dorsal and
        ventral points with normal orthogonal to the line and in the x-y plane.
        """

        def get_ventral_periclinal_points():
            """ Calculate the ventral to periclinal points """

            if thickness_idx == -1 and self.stoma_cfg.wall_thickness.is_thicken_vp_on:
                # thicken cell wall between ventral and periclinal walls

                vp_ellipse = el.SuperEllipse3D( pt_v, pt_p, pt_mid, r=1.5 )  # r=2 => ellipse

                _pts = vp_ellipse.calculate_equi_t_spaced_pts( self.mesh_cfg.num_pts_on_qtr_circumference )

                # TODO: Try this alternative parameterisation of the inner ventral-periclinal curve to get thickening...
                #   (x/a)**2 + (1 + exp(-c * x**2))/2 * (y/b)**2 = 1
                # where 'c' dictates the distance from the ellipse: c=0 => ellipse, c=1 => max deflection.
                #
                # gnuplot:
                # >> a=2; b=1; ai=1.5; bi=0.5
                # >> p [0:2][0:2] for [ al = 0:10 ] bi * sqrt( (1+exp(-0.1*al*x**2))/2 * (1.0 - (x/ai)**2) ) t "".al, \
                #     b * sqrt( 1.0 - (x/a)**2  )
                #
                # The curve is similar to the SuperEllipse but lies closer to the outer ellipse at its ends
                #
                # Note: the super ellipse is: (x/a)**r + (y/b)**r = 1
            else:
                vp_ellipse = el.Ellipse3D( pt_v, pt_p, pt_mid )

                _pts = vp_ellipse.calculate_equi_arc_spaced_pts( self.mesh_cfg.num_pts_on_qtr_circumference )

            return _pts

        def get_periclinal_dorsal_points():
            """ Calculate the periclinal to dorsal points """
            pd_ellipse = el.Ellipse3D( pt_p, pt_d, pt_mid )
            _pts = pd_ellipse.calculate_equi_arc_spaced_pts( self.mesh_cfg.num_pts_on_qtr_circumference )
            return _pts

        # circumferential index of the periclinal points
        periclinal_idx = self.mesh_cfg.num_pts_on_qtr_circumference - 1

        for slice_idx in range( self.mesh_cfg.num_slices ):
            slice_pts = np.ndarray( shape=( self.mesh_cfg.num_pts_through_wall,
                                            self.mesh_cfg.num_pts_on_semi_circumference ),
                                    dtype=p.Point )

            for thickness_idx in ( 0, -1 ):
                # 0 - the outside, -1 - the inside

                # get ventral, periclinal and dorsal points
                pt_v = self.mesh.nodes_grid[ thickness_idx, 0, slice_idx ]
                pt_p = self.mesh.nodes_grid[ thickness_idx, periclinal_idx, slice_idx ]
                pt_d = self.mesh.nodes_grid[ thickness_idx, -1, slice_idx ]

                pt_mid = p.Point( pt_p.x, pt_p.y, 0.0 )

                vp_pts = get_ventral_periclinal_points()
                pd_pts = get_periclinal_dorsal_points()

                # join the circumferential points
                #
                slice_pts[ thickness_idx, : ] = vp_pts[ :-1 ] + pd_pts[ : ]

            # join the points through the wall thickness
            #
            for circumf_idx in range( self.mesh_cfg.num_pts_on_semi_circumference ):
                pt_o = slice_pts[  0, circumf_idx ]
                pt_i = slice_pts[ -1, circumf_idx ]

                line_oi = line.Line( pt_o, pt_i )

                slice_pts[ :, circumf_idx ] = [ line_oi.get_point( t_param )
                                                for t_param in np.linspace( 0.0, 1.0,
                                                                            num=self.mesh_cfg.num_pts_through_wall ) ]

            # add the points to the mesh
            #
            for thickness_idx in range( self.mesh_cfg.num_pts_through_wall ):
                self._add_points_along_circumference( thickness_idx, slice_idx, slice_pts[ thickness_idx, : ] )

        return

    def _add_tip_wall_points( self ):
        """
        The tip wall will be hexahedra with wedges at the centre

        :return:
        """

        def calculate_inner_radial_points():
            """ Calculate points on the inside of the tip wall """

            _pts = list( )

            for _circumf_idx in range( self.mesh_cfg.num_pts_on_semi_circumference ):
                pt_inner_on_wall = tip_wall_nodes[ 0, _circumf_idx, -1 ]
                pt_centre = tip_wall_nodes[ -1, _circumf_idx, -1 ]

                # Use a straight line between the point on the wall and the tip wall centre to
                # space out the points and then drop them onto the wall by setting the y-coordinate
                l_wall_to_centre = line.Line( pt_inner_on_wall, pt_centre )

                inner_pts = list( )
                for _t in np.linspace( 0.0, 1.0, num_radial_pts ):
                    inner_pt = l_wall_to_centre.get_point( _t )
                    inner_pt.y = pt_centre.y
                    inner_pts.append( inner_pt )

                _pts += inner_pts

            return _pts

        def get_interpolator(_pts):
            """ Create bivariate spline interpolator """

            x = [ pt.x for pt in _pts ]
            y = [ pt.y for pt in _pts ]
            z = [ pt.z for pt in _pts ]

            _interp_f = interpolate.SmoothBivariateSpline( x, z, y, kx=4, ky=4 )

            return _interp_f

        # tip_wall_nodes
        # t - radial: 0 (neighbouring the wall)  -->  -1 (the centre)
        # c - circumference: 0 (neighbouring the ventral wall)  -->  -1 (dorsal)
        # s - slice: 0 (shared wall/external)  -->  -1 (internal / pressure surface)
        #
        tip_wall_nodes = self.mesh.tip_nodes_grid

        num_radial_pts = tip_wall_nodes.shape[ 0 ]
        num_pts_through_tip_wall = tip_wall_nodes.shape[ 2 ]

        # copy the inner ventral points into the tip wall nodes grid
        tip_wall_nodes[ 0, :, : ] = self.mesh.nodes_grid[ -1, :, 0:num_pts_through_tip_wall ]

        # set the centre points (through the wall thickness) - repeated in the array
        #
        mid_pts = [ p.Point( self.mesh_cfg.stoma_cfg.mid_ellipse.semi_x_axis, y, 0 )
                    for y in np.linspace( 0.0, self.mesh_cfg.stoma_cfg.wall_thickness.polar,
                                          num=num_pts_through_tip_wall ) ]

        for slice_idx, pt in enumerate( mid_pts ):
            tip_wall_nodes[ -1, :, slice_idx ] = pt

        # add radial points on inside of tip wall using an ellipse
        #
        tip_wall_surface_pts = calculate_inner_radial_points()

        # get an interpolator for the inside of the tip wall
        interp_f = get_interpolator( tip_wall_surface_pts )

        # add points for the tip wall
        #
        for circumf_idx in range( self.mesh_cfg.num_pts_on_semi_circumference ):
            # firstly, add points on the outside of the tip wall (along a line)

            pt_inner = tip_wall_nodes[ 0, circumf_idx, 0 ]
            pt_centre = tip_wall_nodes[ -1, circumf_idx, 0 ]

            l_ic = line.Line( pt_inner, pt_centre )

            t_values = np.linspace( 0.0, 1.0, num_radial_pts )[ 1:-1 ]

            for en_t_idx, t_value in enumerate( t_values ):
                radial_idx = en_t_idx + 1

                pt_tip_wall_outside = l_ic.get_point( t_value )
                # pt_tip_wall_outside = p.from_triad( pt_tip_wall_outside.xyz )

                tip_wall_nodes[ radial_idx, circumf_idx, 0 ] = pt_tip_wall_outside

                # join the outside of the tip wall to the inside

                # calculate the internal tip wall point by interpolation...
                inside_pt_y = interp_f( pt_tip_wall_outside.x, pt_tip_wall_outside.z )

                pt_tip_wall_inside = p.Point( pt_tip_wall_outside.x,
                                              max( self.stoma_cfg.wall_thickness.polar, inside_pt_y[ 0 ] ),
                                              pt_tip_wall_outside.z )

                tip_wall_nodes[ radial_idx, circumf_idx, -1 ] = pt_tip_wall_inside

                # ...and join it to the outside
                l_io = line.Line( pt_tip_wall_outside, pt_tip_wall_inside )

                for en_s_idx, s_value in enumerate( np.linspace( 0.0, 1.0, num_pts_through_tip_wall )[ 1:-1 ] ):
                    slice_idx = en_s_idx + 1

                    pt_in_tip_wall = l_io.get_point( s_value )
                    pt_in_tip_wall = p.from_triad( pt_in_tip_wall.xyz )

                    tip_wall_nodes[ radial_idx, circumf_idx, slice_idx ] = pt_in_tip_wall

        # Overwrite the last tip wall node (that neighbours the GC wall) in case the dorsal wall pt juts into the wall
        for circumf_idx in range( self.mesh_cfg.num_pts_on_semi_circumference ):
            for slice_idx in range( 1, num_pts_through_tip_wall ):
                # move the point that neighbours the wall
                new_pt_1 = sum( tip_wall_nodes[ idx, circumf_idx, slice_idx ] for idx in (0, 2) ) / 2
                tip_wall_nodes[ 1, circumf_idx, slice_idx ] = new_pt_1

                # move the next point too - not always a problem but it can when the stoma is closed
                new_pt_2 = sum( tip_wall_nodes[ idx, circumf_idx, slice_idx ] for idx in (1, 3) ) / 2
                tip_wall_nodes[ 2, circumf_idx, slice_idx ] = new_pt_2

        # add the nodes to the mesh (this will set the node id)
        for nodes_array in (tip_wall_nodes[ 1:, :, : ],):
            for idx_pt_pair in np.ndenumerate( nodes_array ):
                self.mesh.add_node( idx_pt_pair[1] )

        return


    def _add_points_along_axis( self, thickness_idx, circumf_idx, the_pts ):
        """ Add points along an axis """

        for slice_idx, pt in enumerate( the_pts ):
            if self.mesh.nodes_grid[ thickness_idx, circumf_idx, slice_idx ] is None:
                self.mesh.nodes_grid[ thickness_idx, circumf_idx, slice_idx ] = pt
                self.mesh.add_node( pt )

        return self.mesh.nodes_grid[ thickness_idx, circumf_idx, : ]


    def _add_points_along_thickness( self, circumf_idx, slice_idx, the_pts ):
        """ Add points through the wall thickness """

        for thickness_idx, pt in enumerate( the_pts ):
            if self.mesh.nodes_grid[ thickness_idx, circumf_idx, slice_idx ] is None:
                self.mesh.nodes_grid[ thickness_idx, circumf_idx, slice_idx ] = pt
                self.mesh.add_node( pt )

        return self.mesh.nodes_grid[ :, circumf_idx, slice_idx ]


    def _add_points_along_circumference( self, thickness_idx, slice_idx, the_pts ):
        """ Add points along the GC circumference """

        for circumf_idx, pt in enumerate( the_pts ):
            if self.mesh.nodes_grid[ thickness_idx, circumf_idx, slice_idx ] is None:
                self.mesh.nodes_grid[ thickness_idx, circumf_idx, slice_idx ] = pt
                self.mesh.add_node( pt )

        return self.mesh.nodes_grid[ thickness_idx, :, slice_idx ]

    def _add_gc_elements( self ):
        """
        Add the wall elements ensuring the nodes are ordered anti-clockwise for a RH system

        :param mesh:
        :return:
        """

        def calculate_fibre_vector( ):
            """ Calculate the fibre vector """

            ids = ele.node_ids
            p1 = sum( self.mesh.nodes_map[ _ ] for _ in ids[ 0:4 ] ) / 4
            p2 = sum( self.mesh.nodes_map[ _ ] for _ in ids[ 4: ] ) / 4

            return p1 - p2

        # index order is thickness, circumference and then slice
        #
        indices = ((0, 0, 0), (1, 0, 0), (1, 0, 1), (0, 0, 1),
                   (0, 1, 0), (1, 1, 0), (1, 1, 1), (0, 1, 1))

        for thickness_idx in range( self.mesh_cfg.num_pts_through_wall - 1 ):
            for circumf_idx in range( self.mesh_cfg.num_pts_on_semi_circumference - 1 ):
                for slice_idx in range( self.mesh_cfg.num_slices - 1 ):
                    nid_indices = [ (ii[ 0 ] + thickness_idx,
                                     ii[ 1 ] + circumf_idx,
                                     ii[ 2 ] + slice_idx) for ii in indices ]

                    nids = [ self.mesh.nodes_grid[ t, c, s ].id for t, c, s in nid_indices ]

                    ele = elem.Hex8Element( nids )
                    self.mesh.add_element( ele )

                    ele.fibre_vector = calculate_fibre_vector( )

        return

    def _add_tip_wall_elements( self ):
        """
        The polar wall is meshed so that it consists entirely of hexahedral elements.
        This means squares are built at the centre and then these fan out to the PM walls

        :param mesh:
        :param num_pts_circumf:
        :param slice0_pts:
        :param slice1_pts:
        :return:
        """

        def create_element( create_wedge, th_idx, cir_idx, sl_idx ):
            """ Create a tip wall element """

            indices6 = ((0, 0, 0), (0, 1, 0), (1, 0, 0),
                        (0, 0, 1), (0, 1, 1), (1, 0, 1))

            indices8 = ((0, 0, 0), (1, 0, 0), (1, 0, 1), (0, 0, 1),
                        (0, 1, 0), (1, 1, 0), (1, 1, 1), (0, 1, 1))

            indices, ele_class = (indices6, elem.Penta6Element) if create_wedge else (indices8, elem.Hex8Element)

            nid_indices = [ (ii[ 0 ] + th_idx, ii[ 1 ] + cir_idx, ii[ 2 ] + sl_idx)
                            for ii in indices ]

            nids = [ self.mesh.tip_nodes_grid[ t, c, s ].id for t, c, s in nid_indices ]

            ele = ele_class( nids )
            ele.fibre_vector = calculate_fibre_vector( ele )

            return ele

        def calculate_fibre_vector( ele ):
            """ Calculate the vector for the fibre """

            if ele.num_nodes( ) == 8:
                ids = ele.node_ids
                p1 = sum( self.mesh.nodes_map[ _ ] for _ in ids[ 0:4 ] ) / 4
                p2 = sum( self.mesh.nodes_map[ _ ] for _ in ids[ 4: ] ) / 4

                return p1 - p2

            if ele.num_nodes( ) == 6:
                local_facets = ele.local_facets( )
                facet1 = ele.global_facet( local_facets[ 3 ] )
                facet2 = ele.global_facet( local_facets[ 2 ] )

                p1 = sum( self.mesh.nodes_map[ nid ] for nid in facet1 ) / len( facet1 )
                p2 = sum( self.mesh.nodes_map[ nid ] for nid in facet2 ) / len( facet2 )

                return p1 - p2

        # central_square_nodes_grid = mesh.central_square_nodes_grid

        num_radial_pts = self.mesh.tip_nodes_grid.shape[ 0 ]

        # add the wedge elements
        # index order is thickness (inner ventral to tip wall centre), circumference (ventral -> dorsal) and then slice
        #
        for thickness_idx in range( self.mesh.tip_nodes_grid.shape[ 0 ] - 1 ):
            for circumf_idx in range( self.mesh.tip_nodes_grid.shape[ 1 ] - 1 ):
                for slice_idx in range( self.mesh.tip_nodes_grid.shape[ 2 ] - 1 ):
                    # wedges at the centre
                    is6 = (thickness_idx == num_radial_pts - 2)

                    ele = create_element( create_wedge=is6, th_idx=thickness_idx, cir_idx=circumf_idx,
                                          sl_idx=slice_idx )

                    self.mesh.add_element( ele )

        return


def build_stoma_mesh( stoma_cfg, perform_mesh_checks=True ):
    """
    Build stoma mesh with elliptical outline (when observed from above)

    :param stoma_cfg: Stoma configuration object
    :type stoma_cfg: sc.StomateConfig
    :param perform_mesh_checks:
    :return: the mesh
    :rtype: simple_mesh.Mesh
    """

    stoma_mesh = StomaMesh( stoma_cfg )
    stoma_mesh.build( perform_mesh_checks )

    return stoma_mesh.mesh

if __name__ == '__main__':
    pass
