""" A 3D mesh representing a stoma """

from __future__ import print_function

import collections
import logging as lg

import numpy as np
from scipy import sparse
from sortedcontainers import SortedDict

from stomasimulator.geom import point as p

import stomasimulator.stomata.stoma_config as sc

from . import elements as elem
from . import mesh_check
from . import mesh_config

def initialise_ids( ):
    """ Reset all stored identifiers so that they start from 1 """

    # Point
    p.Point.point_id = 1

    elem.BlockElement.element_id = 1


class Mesh( object ):
    """ Represent the 3D mesh of a stoma """

    def __init__( self, stoma_cfg=None ):
        """
        Initialise Mesh
        :param stoma_cfg:
        :type stoma_cfg: sc.StomaConfig
        """

        initialise_ids()

        self.nodes_map = SortedDict( )
        self.edges = list( )
        self.elements = [ ]
        self.bcs = { 'x': [ ], 'y': [ ], 'z': [ ] }

        # generate dicts for these as the key is the guard cell id (0 or 1) and the value is a list of facets
        self.pressure_facets = [ list( ) for _ in range( stoma_cfg.num_gcs ) ]
        self._dorsal_facets = { k: list( ) for k in range( stoma_cfg.num_gcs ) }
        self._ventral_facets = { k: list( ) for k in range( stoma_cfg.num_gcs ) }

        self._epidermis_facets = list( )

        # helper dictionary to store points on the outside and inside
        self._special = dict( node_ids=[ ], pressure_node_ids=[ ], dorsal_node_ids=[ ], ventral_node_ids=[ ],
                              epidermis=[ ] )

        self.mesh_config = mesh_config.MeshConfig( stoma_cfg )

        # structured grids for each part of the mesh
        self._nodes_grids = { k: None for k in ('tube', 'tip', 'epidermis') }

    @property
    def dorsal_wall_facets( self ):
        """ The facets on the dorsal wall """
        return self._dorsal_facets

    @property
    def ventral_wall_facets( self ):
        """ The facets on the vetral wall """
        return self._ventral_facets

    @property
    def epidermis_facets( self ):
        """ The facets on the epidermal wall """
        return self._epidermis_facets

    @property
    def nodes_grid( self ):
        """
        Structured grid of nodes in the tubular part of the stoma (not the tip wall)

        index 0: thickness (from outside to inside)
        index 1: circumference (from ventral to dorsal)
        index 2: slice (from tip to middle)

        :return 3D array
        :rtype np.ndarray
        """

        key = 'tube'

        # structured grid of nodes in the tubular part of the stoma (not the tip wall)
        if self._nodes_grids[ key ] is None:
            mc = self.mesh_config

            shape = (mc.num_pts_through_wall,
                     mc.num_pts_on_semi_circumference,
                     mc.num_slices)

            self._nodes_grids[ key ] = np.ndarray( shape, dtype=p.Point )

        return self._nodes_grids[ key ]

    @property
    def tip_nodes_grid( self ):
        """
        Get the grid of points for the tip wall (initialise it first if necessary)

        index 0: points along the radius of the tip wall (0 is the centre)
        index 1: matches the circumference (from z=0 to z>0)
        index 2: thickness of the tip wall (from y=0 -> +ve, i.e. from outside to inside)

        :return 3D array
        :rtype np.ndarray
        """

        key = 'tip'

        if self.nodes_grid is None:
            raise LookupError( 'The points for the tube must be calculated before the tip wall!' )

        # structured grid of nodes in the tubular part of the stoma (not the tip wall)
        if self._nodes_grids[ key ] is None:
            # get the inner ventral points
            inner_ventral_pts = self.nodes_grid[ -1, 0, : ]

            num_pts_through_tip_wall = self.mesh_config.num_pts_through_tip_wall( inner_ventral_pts )

            shape = ( self.mesh_config.num_pts_tip_wall_radius,
                      self.mesh_config.num_pts_on_semi_circumference,
                      num_pts_through_tip_wall )

            self._nodes_grids[ key ] = np.ndarray( shape, dtype=p.Point )

        return self._nodes_grids[ key ]

    @property
    def epidermal_nodes_grid( self ):
        """
        Structured grid of nodes for the vertical and separate epidermal wall

        index 0: thickness (from neighbouring the GC's dorsal wall outwards)
        index 1: circumference (from z=0 to z>0)
        index 2: slice (from tip to middle)

        :return: 3D array
        :rtype np.ndarray
        """

        key = 'epidermis'

        # structured grid of nodes on the (separate) epidermal wall
        if self._nodes_grids[ key ] is None:
            shape = ( 2,
                      self.mesh_config.num_pts_on_epidermis_vert,
                      2 * self.mesh_config.num_slices)

            self._nodes_grids[ key ] = np.ndarray( shape, dtype=p.Point )

        return self._nodes_grids[ key ]

    def add_node( self, node ):
        """ Add a node to the mesh """

        assert isinstance( node, p.Point ), 'The node is not a Point object!'

        # make sure the node's id is set
        if not node.is_id_set( ):
            node.allocate_id( )

        if node.id in self.nodes_map:
            if self.nodes_map[ node.id ] == node:
                return
            raise ValueError( 'The node id {} already exists in the mesh for a different node!'.format( node.id ) )

        self.nodes_map[ node.id ] = node

        return

    def add_nodes( self, nodes ):
        """
        Add a list of nodes
        :param nodes: list of nodes
        :return:
        """

        if nodes is None:
            return

        for node in nodes:
            self.add_node( node )

    def add_plane_reflection( self, exclude_pt_fn, new_pt_fn, refl_plane ):
        """ Perform a reflection to build the rest of the mesh """

        def create_points():
            """ Create the points """

            # store the new points locally
            new_pts = [ ]

            # create the new points - do not duplicate points (use the exclude_pt_fn function to check)
            for pt in self.nodes_map.values():
                if exclude_pt_fn( pt ):
                    id_map[ pt.id ] = pt.id
                    continue

                new_pt = new_pt_fn( pt )
                new_pts.append( new_pt )

                id_map[ pt.id ] = new_pt.id

                # if the reference point is special then track the reflected one too
                if pt.id in self._special[ 'node_ids' ]:
                    for key in ('pressure_node_ids', 'dorsal_node_ids', 'ventral_node_ids'):
                        if pt.id in self._special[ key ]:
                            self.add_special_node_id( key, new_pt.id )
                            break

            # add the nodes to the mesh
            self.add_nodes( new_pts )

        def calc_fibre_vector( fv ):
            """ Calculate the fibre vector of the new element from the reference fibre vector """

            # don't bother with elements without fibres
            if fv is None:
                return None

            new_fv = p.Point( -fv.x if refl_plane == 'yz' else fv.x,
                              -fv.y if refl_plane == 'xz' else fv.y,
                              -fv.z if refl_plane == 'xy' else fv.z )

            # this negation orients the vector in the same sense
            if refl_plane == 'xy':
                new_fv = -new_fv

            return new_fv

        def create_elements():
            """ Create the elements """

            new_eles = [ ]

            for ele in self.elements:
                # be careful with the node id order
                # - a reflection changes the order and inverts the element...
                node_ids = [ id_map[ node_id ] for node_id in ele.node_ids ]
                num_nids = len( node_ids )
                # ...it's easily handled however by doing this
                node_ids = node_ids[ num_nids // 2: ] + node_ids[ 0:num_nids // 2 ]

                new_ele = elem.BlockElement.build_element_from_node_ids( node_ids )

                new_eles.append( new_ele )

                new_ele.fibre_vector = calc_fibre_vector( ele.fibre_vector )

            self.add_elements( new_eles )

        # map the existing ids to the new ids - used when constructing elements
        id_map = { }

        create_points()
        create_elements()

        return

    def add_special_node_id( self, key, node_id ):
        """ Add a special node """

        self._special[ 'node_ids' ].append( node_id )
        self._special[ key ].append( node_id )

        return

    def add_xy_reflection( self ):
        """ Perform a reflection in the x-y plane"""

        exclude_point = lambda pt: pt.z == 0.0
        new_point = lambda pt: p.Point( pt.x, pt.y, -pt.z, allocate_id=True )

        self.add_plane_reflection( exclude_point, new_point, 'xy' )

        return

    def add_yz_reflection( self ):
        """ Perform a reflection in the y-z plane """

        exclude_point = lambda pt: pt.x == 0.0
        new_point = lambda pt: p.Point( -pt.x, pt.y, pt.z, allocate_id=True )

        self.add_plane_reflection( exclude_point, new_point, 'yz' )

    def add_xz_reflection( self ):
        """ Perform a reflection in the x-z plane """

        exclude_point = lambda pt: pt.y == 0.0
        new_point = lambda pt: p.Point( pt.x, -pt.y, pt.z, allocate_id=True )

        self.add_plane_reflection( exclude_point, new_point, 'xz' )

    def add_edge( self, nid_1, nid_2 ):
        """
        Add an edge to the mesh
        :param nid_1: node id #1
        :param nid_2: node id #2
        :return:
        """

        if nid_1 not in self.nodes_map and nid_2 not in self.nodes_map:
            raise ValueError( 'Tried to add an unknown node id to an edge' )

        self.edges.append( (nid_1, nid_2) )

    def get_adjacency_matrix( self ):
        """
        Get the adjacency matrix for the mesh
        :return: the sparse matrix
        :rtype: np.dok_matrix
        """

        num_nodes = len( self.nodes_map )

        # 1-based - like the point ids
        am = sparse.dok_matrix( (num_nodes + 1, num_nodes + 1), dtype=np.bool )

        for ele in self.elements:
            edges = ele.get_edges( )

            for edge in edges:
                am[ edge[ 0 ], edge[ 1 ] ] = True
                am[ edge[ 1 ], edge[ 0 ] ] = True

        return am

    def add_element( self, element ):
        """ Add an element to the mesh
        :param element:
        :return:
        """

        self.elements.append( element )

        lg.info( 'add_element: %s element added, node ids are %s', element.element_type( ), element.node_ids )

    def add_elements( self, element_list ):
        """ Add a list of elements to the mesh
        :param element_list:
        :return:
        """

        for ele in element_list:
            self.add_element( ele )

    def add_bc_node( self, axis, node_id ):
        """ Set a node to be a boundary condition
        :param axis: x, y or z
        :param node_id:
        :return:
        """
        self.bcs[ axis ].append( node_id )

    def add_pressure_facet( self, gc_index, facet ):
        """ Add a pressure facet (a facet on the inside of a guard cell)
        :param gc_index: 0 or 1
        :param facet: list of node ids
        :return:
        """
        self.pressure_facets[ gc_index ].append( facet )

    def add_external_facet( self, facet ):
        """ Add a facet on the outside of the guard cell to the dorsal/ventral facet list
        :param facet:
        :return:
        """
        # TODO use the node map

        # Is it dorsal or ventral?

        mid_ellipse = self.mesh_config.stoma_cfg.mid_ellipse

        semi_maj = mid_ellipse.semi_x_axis
        semi_min = mid_ellipse.semi_y_axis

        # get the points for the facet
        points = [ self.nodes_map[ nid ] for nid in facet ]

        centroid = sum( points ) / len( points )

        # use the axial ellipse to determine if the facet is dorsal or ventral
        is_ventral_facet = (centroid.x / semi_maj) ** 2 + (centroid.y / semi_min) ** 2 < 1.0

        # GC #0 has y>0 and GC #1 has y<0
        gc_index = 0 if centroid.y > 0 else 1

        if is_ventral_facet:
            self._ventral_facets[ gc_index ].append( facet )
        else:
            self._dorsal_facets[ gc_index ].append( facet )

        return

    def add_epidermis_facet( self, facet ):
        """ Add a facet to the epidermal wall
        :param facet:
        :return:
        """
        self._epidermis_facets.append( facet )

    def get_elements_by_type( self ):
        """ Get elements in the mesh
        :return: dict with keys equal to the element types, e.g. hex8
        """

        result = { }

        for ele in self.elements:
            if ele.element_type( ) not in result:
                result[ ele.element_type( ) ] = [ ]

            result[ ele.element_type( ) ].append( ele )

        return result

    def set_bcs( self ):
        """ Add the boundary conditions (for symmetry) """

        # floating point tolerance
        fp_tolerance = 1e-5

        for node in self.nodes_map.values():
            for coord_label_pair in ((node.x, 'x'), (node.y, 'y'), (node.z, 'z')):
                if abs( coord_label_pair[ 0 ] ) < fp_tolerance:
                    self.add_bc_node( coord_label_pair[ 1 ], node.id )

        return

    def calculate_all_surface_facets( self ):
        """ Calculate information on facets that lie on the surfaces
        When there is 1 guard cell the external tip wall is ignored and when there are 2 GCs it is still ignored."""

        def is_polar_wall_facet( _facet ):
            """
            Facet is on the (external) polar wall if all of its points have x, y or z equal to 0.
            Only happens when mesh has not been reflected, i.e. facet should be part of a wall
            """

            _sum = [ 0.0, 0.0, 0.0 ]
            for _node_id in _facet:
                _sum = [ _sum[ i ] + abs( self.nodes_map[ _node_id ].xyz[ i ] ) for i in range( 3 ) ]

            return abs( min( _sum ) ) < 1e-5

        # map of facets with sorted node ids to facet with the node ids in their original order
        # e.g if facet is (1,3,2) then key-value pair is (1,2,3) -> (1,3,2)
        dict_sortedfacet_facet = dict( )

        # list of facets - each facet has sorted node ids
        # there are duplicates but this is desired!
        #
        all_sortedfacets = [ ]

        for ele in self.elements:
            facets = ele.facets( )
            for facet in facets:
                sortedfacet = tuple( sorted( facet ) )

                # the value will be overwritten when a facet is shared but this occurs when
                # it is in the wall and we are not interested in these ones
                dict_sortedfacet_facet[ sortedfacet ] = facet

                all_sortedfacets.append( sortedfacet )

        # count occurrences of each facet
        sortedfacet_count = collections.Counter( all_sortedfacets )

        # a list of surface facets - the node ids are sorted
        surface_sortedfacets = [ ]

        for candidate_sortedfacet, count in sortedfacet_count.items( ):
            # surface facets have a count of 1 - internal facets are shared so count > 1
            if count == 1 and not is_polar_wall_facet( candidate_sortedfacet ):
                surface_sortedfacets.append( candidate_sortedfacet )

        # remove entries that aren't required
        dict_sortedfacet_facet = { k: dict_sortedfacet_facet[ k ] for k in surface_sortedfacets }

        # create a dictionary of node ids to a list of associated facets
        # e.g. 1 -> [ (1,2,3), (1,4,10) ]
        node_sfacet_dict = dict( )

        for s_sfacet in surface_sortedfacets:
            for node_id in s_sfacet:
                if node_id not in node_sfacet_dict:
                    node_sfacet_dict[ node_id ] = [ ]

                node_sfacet_dict[ node_id ].append( s_sfacet )

        return surface_sortedfacets, dict_sortedfacet_facet, node_sfacet_dict

    @staticmethod
    def calculate_attached_facets( start_node_id, node_sfacet_dict ):
        """ Calculate the facets attached to a node """

        # the set of candidate nodes - will spread from start point like a 'front'
        candidate_nodes = { start_node_id }

        visited_nodes = set( )
        attached_facets = set( )

        while candidate_nodes:
            candidate_node_id = candidate_nodes.pop( )

            if candidate_node_id in visited_nodes:
                continue

            if candidate_node_id in node_sfacet_dict:
                for facet in node_sfacet_dict[ candidate_node_id ]:
                    attached_facets.add( facet )

                    for node_id in facet:
                        if node_id != candidate_node_id and node_id not in visited_nodes:
                            candidate_nodes.add( node_id )

            visited_nodes.add( candidate_node_id )

        return attached_facets

    def set_surface_facets( self ):
        """ Get the surface facets (sorted node ids) and the map of sorted facets -> original facets
        Use 's' to denote sorted """

        surface_facets = self.calculate_all_surface_facets( )

        sfacets_facet_dict, node_sfacet_dict = surface_facets[1], surface_facets[2]

        # find pressure facets (inside each guard cell)
        for idx, node_id in enumerate( self._special[ 'pressure_node_ids' ] ):
            pressure_facets = self.calculate_attached_facets( node_id, node_sfacet_dict )

            # add the pressure facets
            for s_sfacet in pressure_facets:
                facet = sfacets_facet_dict[ s_sfacet ]
                self.add_pressure_facet( idx, facet )

        # find external wall facets by starting at the first dorsal node id -- will then get all external facets
        init_dorsal_nid = self._special[ 'dorsal_node_ids' ][ 0 ]
        external_facets = self.calculate_attached_facets( init_dorsal_nid, node_sfacet_dict )

        # add the facets
        for s_sfacet in external_facets:
            facet = sfacets_facet_dict[ s_sfacet ]
            self.add_external_facet( facet )

        # find ventral wall facets
        # DOES the above find all the facets on the external walls?

        if len( self._special[ 'epidermis' ] ) > 0:
            init_epidermis_nid = self._special[ 'epidermis' ][ 0 ]
            epidermis_facets = self.calculate_attached_facets( init_epidermis_nid, node_sfacet_dict )

            # add the facets
            for s_sfacet in epidermis_facets:
                facet = sfacets_facet_dict[ s_sfacet ]
                self.add_epidermis_facet( facet )

        return

    def __repr__( self ):
        return 'Mesh: {} nodes, {} elements'.format( len( self.nodes_map ), len( self.elements ) )

    def perform_simple_check(self):
        """
        Check the mesh using the Jacobian ratio
        :return:
        """
        j_ratios = [ ele.jacobian_ratio( self.nodes_map ) for ele in self.elements ]

        max_jratio = max( j_ratios )

        num_gt_5 = sum( (1 for _ in j_ratios if _ > 5) )
        num_lt_0 = sum( (1 for _ in j_ratios if _ < 0) )

        print( '--> Mesh statistics:' )
        print( '--> {} elements'.format( len( j_ratios ) ) )
        print( '--> max Jacobian ratio : {}'.format( max_jratio ) )
        print( '--> number > 5         : {}'.format( num_gt_5 ) )
        print( '--> number < 0         : {}'.format( num_lt_0 ) )

    def create_vtk_checker( self ):
        """
        Converts the mesh created by the internal mesh creator to a mesh checker

        :param a_mesh: Instance of a Mesh object (from simple_mesh.py)
        :type a_mesh: simple_mesh.SimpleMesh

        :return:
        """

        points = self.nodes_map.values( )
        elements = [ ele.node_ids for ele in self.elements ]

        return mesh_check.VTKMeshChecker.create_from_points_and_elements( points=points,
                                                                          elements=elements )


def build_test_mesh( ):
    """ Build a tests mesh """

    mesh = Mesh( )

    # create the nodes
    pts = list( )
    pts.append( p.Point( 0, 0, 0 ) )
    pts.append( p.Point( 2, 0, 0 ) )
    pts.append( p.Point( 2, 3, 0 ) )
    pts.append( p.Point( 0, 3, 0 ) )
    pts.append( p.Point( 0, 0, 1 ) )
    pts.append( p.Point( 2, 0, 1 ) )
    pts.append( p.Point( 2, 3, 1 ) )
    pts.append( p.Point( 0, 3, 1 ) )

    pts.append( p.Point( 3, 0, 0 ) )
    pts.append( p.Point( 3, 3, 0 ) )
    pts.append( p.Point( 3, 0, 1 ) )
    pts.append( p.Point( 3, 3, 1 ) )

    p.Point.point_id = 1
    for pt in pts:
        pt.allocate_id( )

    mesh.add_nodes( pts )

    mesh.add_element( elem.Hex8Element( pts[ 0:8 ] ) )
    mesh.add_element( elem.Hex8Element( [ pts[ nid - 1 ] for nid in (2, 9, 10, 3, 6, 11, 12, 7) ] ) )

    # add the BCs
    #
    for axis in ('x', 'y', 'z'):
        for bc_node in (1, 2, 3, 4, 9, 10):
            mesh.add_bc_node( axis, bc_node )

    # Add the pressure facets
    #
    gc_id = 1
    mesh.add_pressure_facet( gc_id, (5, 6, 7, 8) )
    mesh.add_pressure_facet( gc_id, (6, 11, 12, 7) )

    # # Add the wall
    # #
    # mesh.add_polar_wall_facet( (1, 4,  3, 2) )
    # mesh.add_polar_wall_facet( (2, 3, 10, 9) )

    return mesh


def main( ):
    """ Main method """
    mesh = build_test_mesh( )
    mesh.perform_simple_check()


if __name__ == "__main__":
    main( )
