""" Process an XPLT file from an FEBio simulation """

from __future__ import print_function

import collections
import math
import os
import sys
import time

import numpy as np

from sortedcontainers import SortedDict

import stomasimulator.geom.extent as ext
import stomasimulator.geom.geom_utils as geom
import stomasimulator.geom.point as point
import stomasimulator.mesh.mesh_check as mesh_ch

from . import ( xplt_calcs as xcalc,
                xplt_classes as xcls,
                xplt_header as hdr,
                xplt_utils as xplt )

OUTPUT_DELIMITER = '#' * 100


class XpltReaderException( Exception ):
    """ An exception for this module """
    pass


class ModelPart( object ):
    """ Represent part of a model according to the XPLT structure """

    def __init__(self, init_dict=None):
        self.label = None
        self.extent = None
        self.node_ids = set()

        """
        The tracked_points dict is:

        gc0-height:     [ Point, Point ]
        gc0-width:      [ Point, Point ]
        gc1-height:     [ Point, Point ]  -  if the 2nd GC is present
        gc1-width:      [ Point, Point ]  -  if the 2nd GC is present
        pore-boundary:  [ Point, ... ]
        pore-length:    [ Point, Point ]
        pore-width:     [ Point, Point ]
        stoma-width:    [ Point, Point ]
        stoma-length:   [ Point, Point ]
        tracked-points: list of Point objects (all of the above including the pore boundary)
        """
        self.tracked_points = SortedDict()

        # List of AttributeCalculator objects
        self.calculators = list()

        if init_dict is not None:
            self.label = init_dict[ 'label' ] if 'label' in init_dict else self.label
            self.extent = init_dict[ 'extent' ] if 'extent' in init_dict else self.extent
            self.node_ids = init_dict[ 'node_ids' ]if 'node_ids' in init_dict else self.node_ids

        if self.label is None:
            self.label = '<unlabelled>'

        # Store the closed surface(s) and the facets that define them (store by surface_id)
        # will map surface id (e.g. 'surface-1') to a dict with keys ( 'facets', 'node_set', 'edge_set' )
        self.closed_surfaces = collections.OrderedDict( )

    def __repr__(self):
        extra = '' if self.extent is None else 'contained by {}'.format( self.extent )
        return "ModelPart: '{}' has {} nodes {}".format( self.label, len( self.node_ids ), extra )

    def set_tracked_point_list(self):
        """ Store the tracked points in one place - it duplicates the points but simplifies access

        :return:
        """

        # Use a dict to remove duplicates
        pt_dict = SortedDict()

        for pts in self.tracked_points.values():
            for pt in pts:
                if pt.id is not None:
                    pt_dict[ pt.id ] = pt

        pt_list = pt_dict.values()

        self.tracked_points[ 'tracked-points' ] = pt_list

    def print_detail(self):
        print( repr(self) )

        self.print_tracked_nodes()
        self.print_surface_details()

    def print_tracked_nodes(self):
        """ Output information about the tracked nodes """

        print( '' )
        print( 'Tracked node pairs/lists:')

        for key, pts in self.tracked_points.items():
            if len( pts ) == 2:
                print( '  {:20s}: {:40s} and {:40s} distance = {:6.3f}'.format(
                    key, pts[0], pts[1], pts[0].distance( pts[1] ) ) )
            else:
                print( '  {:20s}: {} nodes'.format( key, len( pts ) ) )

    def print_surface_details(self):
        """ Output information about the stored surfaces """

        print( '' )
        print( 'Surfaces:' )

        for surface_id, surface_data in self.closed_surfaces.items():
            if len( surface_data ) == 0:
                continue

            print( '  {:20s}: {} nodes, {} edges and {} facets'.format( surface_id,
                                                                        len( surface_data[ 'node_set' ] ),
                                                                        len( surface_data[ 'edge_set' ] ),
                                                                        len( surface_data[ 'facets'   ] ) ) )


class ModelPartStore(object):
    def __init__(self):
        self._model_parts = []
        self._stoma = None
        self._wall = None

    def print_detail(self):
        for mp in self._model_parts:
            mp.print_detail( )

    @property
    def model_parts( self ):
        """ Get the model parts

        :return: the model parts
        :rtype: Tuple[ModelPart]
        """
        return tuple( self._model_parts )

    @property
    def stoma(self):
        return self._stoma

    @property
    def wall(self):
        return self._wall

    def add_model_part( self, model_part ):
        """ Parts of the model are added during the processing of the domain section.
        The model parts are updated while processing the surfaces.

        :param model_part:
        :type: ModelPart
        :return:
        """

        self._model_parts.append( model_part )

    def set_tracked_point_list( self ):
        """ Update the ModelPart objects

        :return:
        """

        for mp in (self.stoma, self.wall):
            if mp is not None:
                mp.set_tracked_point_list( )

    def update_model_parts( self ):
        """ Store the model parts - assumes there are 1 or 2.

        If a model part has 1 or 2 closed surfaces it is assumed to be the stoma,
        otherwise it is assumed to be the wall.

        This *must* be performed after the surfaces have been processed,
        i.e. do not merge it into 'add_model_part'
        """

        if len( self._model_parts ) not in (1, 2):
            raise XpltReaderException( 'Expected 1 or 2 model parts but got {}!'.format( len( self._model_parts ) ) )

        # assume the model part is the guard cell/stoma when there is only one model part
        if len( self._model_parts ) == 1:
            mp = self._model_parts[ 0 ]
            mp.label = 'stoma'

            if len( mp.closed_surfaces ) == 0:
                # single GC without a tip wall will not have a closed surface so add a 'proxy'
                mp.closed_surfaces['surface-1'] = {}

            self._stoma = mp
        else:
            for mp in self._model_parts:
                num_closed_surfaces = len( mp.closed_surfaces )

                if num_closed_surfaces == 0:
                    mp.label = 'wall'
                    self._wall = mp
                else:
                    if num_closed_surfaces > 2:
                        raise XpltReaderException( 'Found {} closed surfaces but only expected 1 or 2'.format(
                            num_closed_surfaces ) )

                    mp.label = 'stoma'
                    self._stoma = mp

    def print_closed_surface_node_set( self ):
        for mp in self._model_parts:
            for surface_id, value in mp.closed_surfaces.items( ):
                facets = value[ 'facets' ]
                node_set = set( )

                for facet in facets:
                    for node_id in facet:
                        node_set.add( node_id )

                print( 'Node set for surface: {}'.format( surface_id ) )
                print( sorted( list( node_set ) ) )


class GeneralStore( object ):
    """ Multi-purpose container for data extracted by the reader """

    def __init__( self ):
        self.num_nodes = None
        self.materials = []
        self.model_dictionary = dict()
        self.pts_t0 = None
        self.elements = []
        self.model_part_store = ModelPartStore()
        self.states = []

    @property
    def num_gcs(self):
        """ Get the number of guard cells (1 or 2)

        :return: number of guard cells (None if the stoma is undefined)
        :rtype: int
        """
        stoma = self.model_part_store.stoma
        return None if stoma is None else len( stoma.closed_surfaces )

    def print_detail(self):
        """ Readable output of the stored data """

        def print_with_underline( value_, ul_char='=' ):
            """ Helper utility to underline some text
            :param value_:
            :param ul_char:
            :return:
            """

            print( value_ )
            print( ul_char[0] * len( value_ ) )

        print( '' )
        print( OUTPUT_DELIMITER )
        print_with_underline( 'GENERAL STORE:' )
        print( '' )

        print( '{:15s} : {}'.format( 'Number of GCs',   self.num_gcs ) )
        print( '{:15s} : {}'.format( 'Number of nodes', self.num_nodes ) )
        print( '{:15s} : {}'.format( 'Materials',       self.materials ) )
        print( '' )

        print_with_underline( 'Model parts', '-' )
        self.model_part_store.print_detail()
        print( '' )

        print_with_underline( 'State data', '-' )
        for key, value in self.model_dictionary.items():
            print( '{:20s} : {}'.format( xplt.get_tag_for_id( key ), value ) )

        print( '{:20s} : From time = {} --> {}'.format( '{} state(s)'.format( len( self.states ) ),
                                                        self.states[0].time,
                                                        self.states[-1].time ) )
        print( '' )
        print( OUTPUT_DELIMITER )


class ResultsFileWriter(object):
    """ Write the results of processing the XPLT file to a text file. """

    def __init__(self, general_store, metrics ):
        """

        :param general_store:
        :type general_store: GeneralStore
        :param metrics:
        """

        self.general_store = general_store
        self.metrics = metrics

    @staticmethod
    def _write_header( file_handle ):
        file_handle.write( '#\n' )
        file_handle.write( '# State information.\n' )
        file_handle.write( '#\n' )

        # TODO write the model_part_store to the file to get a record of the tracked points
        # TODO write details about each AttributeCalculator

    @staticmethod
    def _write_column_headers( file_handle, column_names ):
        """
        Write the column numbers on one line and the column names on the next.

        :param file_handle:
        :param column_names:
        :return:
        """

        # only want comment character before column numbers
        file_handle.write( '# ' )

        # write a 1-based list of column numbers on one line followed by a list of column names
        for coll in ( range( 1, 1 + len( column_names ) ), column_names ):
            for item in coll:
                file_handle.write( '{}  '.format( item ) )

            file_handle.write( '\n' )

    @staticmethod
    def _write_column_data( file_handle, state, column_names,
                            write_time=True, new_line=True, check_names_exist=False ):

        if write_time:
            file_handle.write( '{:f} '.format( state.time ) )

        for column_name in column_names:
            if check_names_exist:
                datum = state.get_attribute( column_name ) if state.has_attribute( column_name ) else float('nan')
            else:
                datum = state.get_attribute( column_name )

            # fmt_string = ' {:12.6f} '
            fmt_string = ' {:12.8f} '
            if isinstance( datum, ( str, int ) ):
                fmt_string = ' {} '

            file_handle.write( fmt_string.format( datum ) )

        if new_line:
            file_handle.write( '\n' )

    def _write_mesh_quality( self, file_handle ):

        file_handle.write( '#\n' )

        for state in self.general_store.states:
            if state.has_attribute( 'mesh-quality' ):
                mesh_quality = state.get_attribute( 'mesh-quality' )

                msg = '# Mesh quality for state {} (t = {})\n'.format( state.ordinal, state.time )
                file_handle.write( msg )

                for metric_quality in mesh_quality:
                    file_handle.write( '# {}\n'.format( metric_quality ) )

                file_handle.write( '#\n' )

    def _write_stoma_data( self, file_handle, dataset_index ):
        """ Write the results for the whole stoma

        :param file_handle:
        :return:
        """

        file_handle.write( '\n\n' )
        file_handle.write( '# dataset {}: For the stoma (column numbers are in brackets)\n'.format( dataset_index ) )

        # fetch state attributes with these keys
        #
        stoma_keys = ( 'pore-length', 'pore-width', 'stoma-length', 'stoma-width', 'pore-area' )

        # guard cell keys are: gc0-height, gc0-width, ... , gc1-height, gc1-width, ...
        #
        gc_key_suffices = ( 'height', 'width', 'ar', 'surface-area', 'volume' )
        gc_keys = []
        for gc_key in ( 'gc0', 'gc1' ):
            gc_keys.append( tuple( '{}-{}'.format( gc_key, key ) for key in gc_key_suffices ) )

        self._write_column_headers( file_handle, ('t',) + stoma_keys + gc_keys[0] + gc_keys[1] )

        # write the data
        #
        for state in self.general_store.states:
            self._write_column_data( file_handle, state, stoma_keys, new_line=False )
            self._write_column_data( file_handle, state, gc_keys[0], write_time=False, check_names_exist=True, new_line=False )
            self._write_column_data( file_handle, state, gc_keys[1], write_time=False, check_names_exist=True )

    def _identify_best_metric(self, best_key, metric_name ):
        states = self.general_store.states

        if len( states ) == 0 or not states[0].has_attribute( metric_name ):
            return

        for state in states:
            state.set_attribute( best_key, 'no' )

        best_state = min( (state.get_attribute( metric_name ), state) for state in states )[ 1 ]
        best_state.set_attribute( best_key, 'BEST-STATE' )

    def _write_metrics(self, file_handler, dataset_index ):

        # do nothing if the optimisation keys are absent
        init_state = self.general_store.states[ 0 ]

        if not init_state.has_attribute( 'optimisation-keys' ):
            return

        file_handler.write( '\n\n' )
        file_handler.write( '# dataset {}: optimisation metric\n'.format( dataset_index ) )

        best_key = 'best'

        keys = init_state.get_attribute( 'optimisation-keys' )
        keys.append( best_key )

        self._write_column_headers( file_handler, ('t',) + tuple( keys ) )

        self._identify_best_metric( best_key, 'metric' )

        # write the data to the file
        #
        for state in self.general_store.states:
            self._write_column_data( file_handler, state, keys, check_names_exist=True )

    def write_to_file(self, file_name ):
        """ Write the results to the file

        :param file_name:
        """

        with open( file_name, 'w' ) as fh:
            self._write_header( fh )
            self._write_mesh_quality( fh )
            self._write_stoma_data( fh, dataset_index=0 )
            self._write_metrics( fh, dataset_index=1 )

        print( 'Written output (area, volume and nodal displacements) to {}'.format( file_name ) )


class XpltFileReader(object):
    """ Read and process the XPLT file """

    def __init__(self, xplt_filename, verbose, metrics ):
        """
        :param xplt_filename:
        :param verbose:
        :param metrics:
        :type metrics: xcalc.XpltReaderMetrics
        """

        self.xplt_filename = xplt_filename
        self.verbose = verbose
        self.metrics = metrics

        self.pore_centre = point.Point( 0, 0, 0 )
        self.general_store = GeneralStore()

        print( OUTPUT_DELIMITER )
        print( '--> XPLT File            : {}'.format( self.xplt_filename ) )
        print( '--> Compute metric       : {}'.format( self.metrics.is_compare_vs_open_stoma_on ) )
        print( '--> Check mesh quality   : {}'.format( self.metrics.is_mesh_calculation_on ) )
        print( '--> Pore centre          : {}'.format( self.pore_centre ) )

    @staticmethod
    def output_raw_values( input_data, offset=0, num=None ):
        """
        Help for debugging - readable output of the binary data
        :param input_data:
        :param offset:
        :param num:
        :return:
        """

        print( OUTPUT_DELIMITER )

        if num is None:
            num = len( input_data ) - offset

        num_col_out = 0
        msg = False

        # Special processing for PLT_FACE and PLT_ELEMENT entries because there are (potentially)
        # thousands of individual entries and the output gets swamped

        num_plt_face, tag_plt_face = 0, hdr.FEBIO_TAGS[ 'PLT_FACE' ]
        num_plt_element, tag_plt_element = 0, hdr.FEBIO_TAGS[ 'PLT_ELEMENT' ]

        tag_plt_state = hdr.FEBIO_TAGS[ 'PLT_STATE' ]

        go_silent = False

        print( 'Raw data...' )
        print( 'Byte-offset Tag                          Tag id      Size       (    size)   Data...' )
        print( '' )

        current_tag_id = None

        problem_tags = ( hdr.FEBIO_TAGS[ 'PLT_SURFACE_NAME' ],
                         hdr.FEBIO_TAGS[ 'PLT_DOM_NAME' ],
                         hdr.FEBIO_TAGS[ 'PLT_NODESET_SECTION' ] )

        xx = offset
        while xx < offset + num - hdr.SZ_DWORD:
            incremented_x = False

            dw = xplt.convert_4bytes_to_dword( input_data, xx )

            if dw in hdr.TAGS:
                current_tag_id = dw

                if dw == tag_plt_face or dw == tag_plt_element:
                    if dw == tag_plt_face:
                        num_plt_face += 1

                    if dw == tag_plt_element:
                        num_plt_element += 1

                    if num_plt_face == 5 or num_plt_element == 5:
                        sys.stdout.write( '\n ...' )
                        go_silent = True
                else:
                    num_plt_face = 0
                    num_plt_element = 0
                    go_silent = False

                if not go_silent:
                    if dw == tag_plt_state:
                        sys.stdout.write( '\n{:s}'.format( OUTPUT_DELIMITER ) )

                    sys.stdout.write( '\n {:10} {:25s} : '.format( xx, xplt.get_tag_for_id( dw ) ) )

                num_col_out = 0
                msg = False

            if not go_silent:
                if num_col_out < 8:
                    sys.stdout.write( ' {:#010x} '.format( dw ) )
                    num_col_out += 1

                    # put the size in brackets
                    if num_col_out == 2:
                        sys.stdout.write( '({:8})  '.format( dw ) )

                        # some names that should be CHAR64 are not - so deal with it
                        if dw % hdr.SZ_DWORD != 0 and current_tag_id in problem_tags:
                            # move past the size
                            xx += hdr.SZ_DWORD

                            # output the byte data - don't convert it to DWORDS - it's not aligned
                            print( ' {}'.format( input_data[ xx: xx + dw ] ), end='' )

                            # move past the byte data
                            xx += dw
                            incremented_x = True
                else:
                    if not msg:
                        sys.stdout.write( '...<truncated>' )
                        msg = True

            if not incremented_x:
                xx += hdr.SZ_DWORD

        print( '' )
        print( '>(end)' if (num + offset == len( input_data )) else '>...' )

        print( OUTPUT_DELIMITER )

    def process_node_section( self, coords ):
        """ Store the node coordinates (converting them to a list of points along the way)

        :param coords: the array of float values
        :type coords: np.ndarray
        """

        xyz = np.frombuffer( coords, dtype=np.float32, offset=0 )

        nn = self.general_store.num_nodes

        # x3 because array stores sets of x,y,z
        #
        if len( xyz ) != nn * 3:
            print( '' )
            print( '* Argh! Unexpected number of nodes' )
            print( '* GEN_STORE[ "PLT_HDR_NODES" ] = {}'.format( nn ) )
            print( '* len( buffer )                = {}'.format( len( xyz ) ) )
            print( '' )

            raise XpltReaderException( 'Encountered an unexpected number of nodes' )

        pts_dict = self.xyz_list_to_dict( xyz )

        # store the coordinates (in xyz-1, xyz-2, ... format)
        self.general_store.pts_t0 = pts_dict

        print( '--> stored coordinates for the {} nodes'.format( nn ) )

    @staticmethod
    def xyz_list_to_dict( flt_list ):
        """ Convert a list of floats to a dict with id -> Point mapping.

        :param flt_list: list of floating point values
        :return: mapping from point id to Point objects
        :rtype: dict(int, Point)
        """

        # start_time = time.time()

        # A node-id to Point dictionary
        pts_dict = dict( )

        len_list = len( flt_list )

        nn = len_list // 3

        if 3 * nn != len_list:
            print( 'process_xyz_list: Error - list length is not a multiple of 3 (length = {})'.format( len_list ) )
            quit( )

        # we rely on this so reset it
        point.Point.point_id = 1

        for x_idx in range( 0, len_list, 3 ):
            xx, yy, zz = flt_list[ x_idx:x_idx + 3 ]
            pt = point.Point( xx, yy, zz, allocate_id=True )
            pts_dict[ pt.id ] = pt

        # print( '--> xyz_list_to_dict: time taken = {}'.format( time.time() - start_time ) )

        return pts_dict

    def process_domain_section( self, domain_section_data ):
        """
        Example data...

        PLT_DOMAIN_SECTION       0x01042000  0x000e3954
        PLT_DOMAIN               0x01042100  0x000e394c
        PLT_DOMAIN_HDR           0x01042101  0x00000024
        PLT_DOM_ELEM_TYPE        0x01042102  0x00000004  0x00000002
        PLT_DOM_MAT_ID           0x01042103  0x00000004  0x00000001
        PLT_DOM_ELEMS            0x01032104  0x00000004  0x0000820a



        PLT_DOM_ELEM_LIST        0x01042200  0x000e3918
        PLT_ELEMENT              0x01042201  0x00000014  0x00000001  0x00001b6c  0x000008e5  0x000009b4  0x000008a4
        PLT_ELEMENT              0x01042201  0x00000014  0x00000002  0x00002679  ...
        ...
        PLT_DOMAIN ...

        PLT_DOMAIN_SECTION        :  0x01042000  0x000385ba (  230842)
        PLT_DOMAIN                :  0x01042100  0x00037b45 (  228165)
        PLT_DOMAIN_HDR            :  0x01042101  0x00000035 (      53)
        PLT_DOM_ELEM_TYPE         :  0x01042102  0x00000004 (       4)   0x00000000
        PLT_DOM_MAT_ID            :  0x01042103  0x00000004 (       4)   0x00000001
        PLT_DOM_ELEMS             :  0x01032104  0x00000004 (       4)   0x00001440
                                     0x01032105  0x00000009  0x00000005  0x74726150  0x04220031  ...

        :param domain_section_data:
        """

        if self.verbose:
            self.output_raw_values( domain_section_data, 0, 1200 )

        # Will be a list of sets of points. Each item in the list will be connected and so define a connected region.
        list_of_node_sets = [ ]

        # Iterate through the PLT_DOMAIN entries (each domain will consist of one type of element and can be a
        # sub-domain of a larger domain)...
        #
        ele_map = { 0: 'HEX8', 1: 'PENTA6', 2: 'TET4', 3: 'QUAD4', 4: 'TRI3', 5: 'TRUSS2' }
        offset = 0
        while True:
            # Find PLT_DOMAIN
            offset = xplt.get_offset_for_tag( domain_section_data, 'PLT_DOMAIN', offset )

            if offset is None:
                break

            domain, offset = self.read_section( domain_section_data, 'PLT_DOMAIN', offset, True )

            domain_data = domain.get_datum_ndarray( )

            offset2 = 0
            dom_hdr, offset2 = self.read_section( domain_data, 'PLT_DOMAIN_HDR', offset2 )
            dom_ele_t, offset2 = self.read_section( domain_data, 'PLT_DOM_ELEM_TYPE', offset2, include_data=True )
            dom_mat_id, offset2 = self.read_section( domain_data, 'PLT_DOM_MAT_ID', offset2, include_data=True )
            dom_elements, offset2 = self.read_section( domain_data, 'PLT_DOM_ELEMS', offset2, include_data=True )

            # Older files do not have the tag
            offset_temp = xplt.get_offset_for_tag( domain_data, 'PLT_DOM_NAME', offset2 )
            if offset_temp is not None:
                dom_name, offset2 = self.read_section( domain_data, 'PLT_DOM_NAME', offset_temp, True )
                print( '--> Domain name is : {}'.format( dom_name.get_datum_str( default='<unnamed>' ) ) )

            dom_elem_list, offset2 = self.read_section( domain_data, 'PLT_DOM_ELEM_LIST', offset2 )

            ele_type_id = dom_ele_t.get_datum_int()
            element_type_str = ele_map.get( ele_type_id ) if ele_type_id in ele_map else '<Unrecognised element type>'

            print("--> found domain of {} {} elements of material {}".format(
                dom_elements.get_datum(), element_type_str, dom_mat_id.get_datum() ) )

            # Get the connectivity set for the points
            node_set = set()
            list_of_node_sets.append( node_set )

            for _ in range( dom_elements.get_datum_int( ) ):
                dom_element, offset2 = self.read_section( domain_data, 'PLT_ELEMENT', offset2, include_data=True )

                el_data = dom_element.get_datum_ndarray( )

                # The element id (stored in el_data[0]) is ignored
                # Add 1 because node ids are 0-based here and 1-based everywhere else
                nids_list = [ node_idx + 1 for node_idx in el_data[1:] ]

                # add the node ids to the set
                node_set |= set( nids_list )

                if self.metrics.is_mesh_calculation_on:
                    self.general_store.elements.append( nids_list )

            if self.verbose:
                print( domain )
                print( dom_hdr )
                print( dom_ele_t )
                print( dom_mat_id )
                print( dom_elements )
                print( dom_elem_list )

        num_node_sets = len( list_of_node_sets )

        print( '--> number of node sets = {}'.format( num_node_sets ) )

        if num_node_sets > 1:
            # Find the node sets that share nodes and merge them to get a list of distinct node sets

            merging = True
            while merging and num_node_sets > 1:
                num_node_sets = len( list_of_node_sets )
                merging = False

                for i in range( num_node_sets ):
                    if list_of_node_sets[ i ] is None:
                        continue

                    for j in range( i + 1, num_node_sets ):
                        if list_of_node_sets[ j ] is None:
                            continue

                        if not list_of_node_sets[ i ].isdisjoint( list_of_node_sets[ j ] ):
                            list_of_node_sets[ i ].update( list_of_node_sets[ j ] )
                            list_of_node_sets[ j ] = None
                            merging = True

                            print( '--> merged node sets {} and {}'.format( i, j ) )

                # Get rid of the blank lists (start at end and pop the 'None's)
                i = num_node_sets - 1
                while i >= 0:
                    if list_of_node_sets[ i ] is None:
                        list_of_node_sets.pop( i )
                    i -= 1

        print( '--> identified {} disjoint node set(s)'.format( len( list_of_node_sets ) ) )

        for node_set in list_of_node_sets:
            an_ext = ext.BoxExtent.from_points( [ self.general_store.pts_t0[ node_id ] for node_id in node_set ] )

            mp = ModelPart( dict( node_ids=node_set, extent=an_ext ) )

            self.general_store.model_part_store.add_model_part( mp )

            print( '--> stored {}'.format( mp ) )

    def process_surface_section( self, surface_section_data ):
        """ Extract the surfaces and store the closed one(s).
        Example data...

        PLT_SURFACE_SECTION      0x01043000  0x000d16e0
        PLT_SURFACE              0x01043100  0x000443b8
        PLT_SURFACE_HDR          0x01043101  0x00000018
        PLT_SURFACE_ID           0x01043102  0x00000004  0x00000001
        PLT_SURFACE_FACES        0x01043103  0x00000004  0x00001b4a
        PLT_FACE_LIST            0x01043200  0x00044390
        PLT_FACE                 0x01043201  0x00000020  0x00000001  0x00000003  0x0000156e  ...
        PLT_FACE                 0x01043201  0x00000020  0x00000002  0x00000003  0x00001824  ...
        ...
        PLT_SURFACE ...

        :param surface_section_data:
        """

        if self.verbose:
            self.output_raw_values( surface_section_data, 0, 1200 )

        # Iterate through the PLT_SURFACE entries...
        #
        offset = 0
        while True:
            # Find PLT_SURFACE
            offset = xplt.get_offset_for_tag( surface_section_data, 'PLT_SURFACE', offset )

            if offset is None:
                break

            p_surface, offset = self.read_section( surface_section_data, 'PLT_SURFACE', offset )
            p_surf_hdr, offset = self.read_section( surface_section_data, 'PLT_SURFACE_HDR', offset )
            p_surf_id, offset = self.read_section( surface_section_data, 'PLT_SURFACE_ID', offset, include_data=True )

            p_surf_faces, offset = self.read_section( surface_section_data,
                                                      'PLT_SURFACE_FACES',
                                                      offset,
                                                      include_data=True )

            # Older files do not have the tag
            offset_temp = xplt.get_offset_for_tag( surface_section_data, 'PLT_SURFACE_NAME', offset )

            if offset_temp is not None:
                p_surf_name, offset = self.read_section( surface_section_data,
                                                         'PLT_SURFACE_NAME',
                                                         offset_temp,
                                                         include_data=True )

                print( '--> Surface name is : {}'.format( p_surf_name.get_datum_str( default='<unnamed>' ) ) )

            p_surf_face_list, offset = self.read_section( surface_section_data, 'PLT_FACE_LIST', offset )

            if self.verbose:
                print( p_surface )
                print( p_surf_hdr )
                print( p_surf_id )
                print( p_surf_faces )
                print( p_surf_face_list )

            surface_id = p_surf_id.get_datum( )
            num_faces = p_surf_faces.get_datum_int( )

            surface = [ ]

            # Read PLT_FACE entries - 10 DWORDS (2 control [ID, num nodes] + 8 data [node ids])
            # e.g.
            #
            # PLT_FACE    SIZE        ID          NUM_NODES   NODE IDs...
            # 0x01043201  0x00000020  0x00000001  0x00000003  0x0000156e  0x00001745  0x00001742  ...
            # 3 nodes     #1          #2          #3          ignore      ignore      ignore
            #
            for _ in range( num_faces ):
                # Get PLT_FACE and its data
                p_face, offset = self.read_section( surface_section_data, 'PLT_FACE', offset, True )

                p_face_data = p_face.get_datum_ndarray( )

                # Index 0 - spec. says to ignore it and it *DOES NOT* match the id in PostView)
                internal_id = p_face_data[ 0 ]

                # index 1 stores the number of nodes
                num_nodes = p_face_data[ 1 ]

                # WARNING Add 1 because the facets appear to store the index rather than the id (which will be off by 1)
                node_ids = tuple( [ x + 1 for x in p_face_data[ 2:2 + num_nodes ] ] )

                if self.verbose:
                    print( 'Internal id {}: Facet has {} nodes with node ids {}'.format(
                        internal_id, num_nodes, node_ids ) )

                surface.append( node_ids )

            # Only store closed surfaces (because we want volume and it's undefined for an open surface)
            chi, node_set, edge_set = geom.calculate_euler_characteristic( surface )

            if chi == 2:
                # The surface is closed => store it
                for mp in self.general_store.model_part_store.model_parts:
                    if not node_set.isdisjoint( mp.node_ids ):
                        surface_key = 'surface-{}'.format( surface_id )

                        mp.closed_surfaces[ surface_key ] =\
                            { 'facets': surface, 'node_set': node_set, 'edge_set': edge_set }

                        print( "--> Stored closed surface {} using key '{}'".format( surface_id, surface_key ) )
            else:
                print( '--> Ignoring surface (id={}) because it is not closed'.format( surface_id ) )

        self.general_store.model_part_store.update_model_parts()

        if self.verbose:
            self.general_store.model_part_store.print_closed_surface_node_set()

    def process_root_section( self, root_section ):
        """
        Read the top-level section and break it down into its constituent parts

        :param root_section:
        :return:
        """

        offset = 0
        root_data = root_section.get_datum( )

        if self.verbose:
            self.output_raw_values( root_data, 0, 600 )
            print( OUTPUT_DELIMITER )

        ##########################
        # Some root header specific processing functions
        #
        def process_header_section_data( hdr_data ):
            """ Shred the header information - get the number of nodes """

            print( '--> Processing root header section...' )

            r_offset = 0

            hdr_version, r_offset = self.read_section( hdr_data, 'PLT_HDR_VERSION', r_offset )
            hdr_nodes, r_offset = self.read_section( hdr_data, 'PLT_HDR_NODES', r_offset )
            hdr_facet, r_offset = self.read_section( hdr_data, 'PLT_HDR_MAX_FACET_NODES', r_offset, True )

            if self.verbose:
                print( hdr_version )
                print( hdr_nodes )
                print( hdr_facet )

            self.general_store.num_nodes = hdr_nodes.get_datum( )

            print( '--> {} nodes'.format( hdr_nodes.get_datum( ) ) )

        def process_dictionary_section_data( d_data ):
            """ Shred the meta-data relating to the model, e.g. what is stored """

            print( '--> Processing dictionary section data...' )

            dict_tags = ( 'PLT_DIC_GLOBAL', 'PLT_DIC_MATERIAL', 'PLT_DIC_NODAL', 'PLT_DIC_DOMAIN', 'PLT_DIC_SURFACE' )

            for dict_tag in dict_tags:
                d_offset = xplt.get_offset_for_tag( d_data, dict_tag, 0 )

                print( '--> is {} section present? {}'.format( dict_tag, 'no ' if d_offset is None else 'yes' ) )
                if d_offset is None:
                    continue

                dict_tag_section, dummy = self.read_section( d_data, dict_tag, d_offset, True )

                # Process the PLT_DIC_ITEMs
                di_offset = 0
                while True:
                    dict_item, di_offset = self.read_section( dict_tag_section.get_datum_ndarray( ),
                                                              'PLT_DIC_ITEM',
                                                              di_offset )

                    if dict_item is None:
                        break

                    dict_item_type, di_offset = self.read_section( dict_tag_section.get_datum_ndarray( ),
                                                                   'PLT_DIC_ITEM_TYPE',
                                                                   di_offset,
                                                                   include_data=True )

                    dict_item_fmt, di_offset = self.read_section( dict_tag_section.get_datum_ndarray( ),
                                                                  'PLT_DIC_ITEM_FMT',
                                                                  di_offset,
                                                                  include_data=True )

                    dict_item_name, di_offset = self.read_section( dict_tag_section.get_datum_ndarray( ),
                                                                   'PLT_DIC_ITEM_NAME',
                                                                   di_offset,
                                                                   include_data=True )

                    dict_item_obj = xcls.DictionaryItem( dict_item_type.get_datum_int(),
                                                         dict_item_fmt.get_datum_int(),
                                                         dict_item_name.get_datum_ndarray() )

                    state_tag_id = hdr.FEBIO_DICT_STATE_MAP.get( dict_tag_section.id )

                    if state_tag_id not in self.general_store.model_dictionary.keys( ):
                        self.general_store.model_dictionary[ state_tag_id ] = [ ]

                    self.general_store.model_dictionary[ state_tag_id ].append( dict_item_obj )

                    print( '--> {}'.format( dict_item_obj ) )

        def process_material_section_data( material_data ):
            """ Shred the data relating to the materials in the model """

            print( '--> Processing material section...' )

            m_offset = 0
            while True:
                material, m_offset = self.read_section( material_data, 'PLT_MATERIAL', m_offset )

                if material is None:
                    break

                mat_id, m_offset = self.read_section( material_data, 'PLT_MAT_ID', m_offset, True )
                mat_name, m_offset = self.read_section( material_data, 'PLT_MAT_NAME', m_offset, True )

                mat_obj = xcls.Material( mat_id.get_datum( ), mat_name.get_datum( ) )

                self.general_store.materials.append( mat_obj )

                print( '--> {}'.format( mat_obj ) )

        def process_geometry_section_data( geom_data ):
            """ Shred the data relating to the geometry """

            print( '--> Processing geometry section...' )

            g_offset = 0

            # Process the node section
            #
            print( '--> nodes...' )

            g_offset = self.read_section( geom_data, 'PLT_NODE_SECTION', g_offset )[ 1 ]
            node_coords, g_offset = self.read_section( geom_data, 'PLT_NODE_COORDS', g_offset, True )

            self.process_node_section( node_coords.get_datum_ndarray( ) )

            # Process domain section
            #
            print( '--> domains...' )

            domain_section, g_offset = self.read_section( geom_data, 'PLT_DOMAIN_SECTION', g_offset, True )

            self.process_domain_section( domain_section.get_datum( ) )

            # Process surface section
            #
            print( '--> surfaces...' )

            # the domain section may have had its size altered (due to the DWORD realignment) so find the next id
            # g_offset = xplt.get_offset_for_tag( geom_data, 'PLT_SURFACE_SECTION', g_offset, self.verbose )

            surface_section, g_offset = self.read_section( geom_data, 'PLT_SURFACE_SECTION', g_offset, True )

            self.process_surface_section( surface_section.get_datum( ) )

            # Find tracked nodes/facets
            #
            print( '--> finding nodes to track...' )
            self.store_nodes_of_interest( )
        #
        # End of root header specific processing functions
        ##########################

        # Process the header
        #
        header_section, offset = self.read_section( root_data, 'PLT_HEADER', offset, True )

        if self.verbose:
            print( header_section )

        process_header_section_data( header_section.get_datum( ) )

        # Process the dictionary
        #
        dict_section, offset = self.read_section( root_data, 'PLT_DICTIONARY', offset, True )

        if self.verbose:
            print( dict_section )

        process_dictionary_section_data( dict_section.get_datum( ) )

        # Process the material(s)
        #
        material_section, offset = self.read_section( root_data, 'PLT_MATERIALS', offset, True )

        if self.verbose:
            print( material_section )

        process_material_section_data( material_section.get_datum( ) )

        # Process the geometry
        #
        geometry_section, offset = self.read_section( root_data, 'PLT_GEOMETRY', offset, True )

        if self.verbose:
            print( geometry_section )

        process_geometry_section_data( geometry_section.get_datum( ) )

    def store_nodes_of_interest( self ):
        """ The nodes of interest are (assuming axis is x, pore/stoma width is y and height is z):
        1. Central points (x=z=0) on both GCs to measure the pore width and the stoma width.
           Points on the cell wall and the PM are stored.
        2. Points that form the boundary of the pore in the equatorial plane (z=0) on the PM ventral wall
        3. Points at the ends to measure the stoma length
        """

        def store_width_points():
            """
            Get the GC points that will be used to measure the pore width and stoma width.
            In both cases we obtain points on the cell wall and the PM.

            Algorithm:
            o  Construct thin cuboid that will contain the points

            NOTE: Assume axis is x!
            :return:
            """

            # build an extent that surrounds the required width points
            main_ext = ext.BoxExtent.from_point( self.pore_centre ).grow( dx=padding, dy=1e6, dz=padding )

            # get the points in the extent
            pts = main_ext.get_contained_points( ( self.general_store.pts_t0.get( _nid ) for _nid in stoma.node_ids ) )

            # Always have 1 guard cell (all points in GC #1 have y >= 0)
            gc_pts = [ _pt for _pt in pts if _pt.y > 0.0 ]
            gc_pts.sort( key=lambda _pt: self.pore_centre.distance( _pt ) )

            pt0_ventral, pt0_dorsal = gc_pts[ 0 ], gc_pts[ -1 ]
            pt1_ventral, pt1_dorsal = self.pore_centre, self.pore_centre

            if self.general_store.num_gcs == 2:
                # Second GC has all of its points in y <= 0
                gc_pts = [_pt for _pt in pts if _pt.y < 0.0]
                gc_pts.sort(key=lambda _pt: self.pore_centre.distance(_pt))

                pt1_ventral, pt1_dorsal = gc_pts[0], gc_pts[-1]

            # store the points
            stoma.tracked_points[ 'stoma-width' ] = ( pt0_dorsal,  pt1_dorsal )
            stoma.tracked_points[ 'pore-width'  ] = ( pt0_ventral, pt1_ventral )
            stoma.tracked_points[ 'gc0-width'   ] = ( pt0_dorsal,  pt0_ventral )

            if self.general_store.num_gcs == 2:
                stoma.tracked_points['gc1-width'] = ( pt1_dorsal, pt1_ventral )

        def store_length_points():
            """
            Get the GC points that will be used to measure the stoma's length (only cell wall points)

            Algorithm:
            o  Construct thin cuboid that will contain the points

            NOTE: Assumes length axis is x!
            :return:
            """

            # Is the padding small enough?

            # build an extent that surrounds the required points
            main_ext = ext.BoxExtent.from_point( self.pore_centre ).grow( dx=1e6, dy=padding, dz=padding)

            # get the points in the extent
            pts = main_ext.get_contained_points( ( self.general_store.pts_t0.get( _nid ) for _nid in stoma.node_ids ) )

            # Order by x coordinate
            pts.sort( key=lambda _pt: _pt.x )

            pos_x_pts = [ _pt for _pt in pts if _pt.x >= 0.0 ]
            neg_x_pts = [ _pt for _pt in pts if _pt.x <= 0.0 ]

            stoma.tracked_points[ 'pore-length'  ] = ( neg_x_pts[-1], pos_x_pts[0] )
            stoma.tracked_points[ 'stoma-length' ] = ( pts[0], pts[-1] )

        def store_height_points():
            """ Find and store the points that track the GC height(s)

            :return:
            """

            # build an extent that contains the required points
            main_ext = ext.BoxExtent.from_point( self.pore_centre ).grow( dx=padding, dy=1e6, dz=1e6 )

            # get the points in the extent
            pts = main_ext.get_contained_points( ( self.general_store.pts_t0.get( _nid ) for _nid in stoma.node_ids ) )

            pts.sort( key=lambda _pt: _pt.z )

            pos_y_pts = [ _pt for _pt in pts if _pt.y >= 0.0 ]
            neg_y_pts = [ _pt for _pt in pts if _pt.y <= 0.0 ]

            stoma.tracked_points[ 'gc0-height' ] = ( pos_y_pts[0], pos_y_pts[-1] )

            if self.general_store.num_gcs == 2 and len( neg_y_pts ) > 1:
                stoma.tracked_points[ 'gc1-height' ] = ( neg_y_pts[0], neg_y_pts[-1] )

        def store_pore_boundary_nodes( ):
            """ Find and store the nodes that are used for pore-area calculations """

            # build an extent that contains the required points
            main_ext = ext.BoxExtent.from_point( self.pore_centre ).grow( dx=1e6, dy=1e6, dz=padding )

            pts = main_ext.get_contained_points( ( self.general_store.pts_t0.get( _nid ) for _nid in stoma.node_ids) )

            # filter using the pore ellipse, i.e. get rid of the points that lie on or outside the ellipse
            semi_maj = self.pore_centre.distance( stoma.tracked_points[ 'pore-length' ][0] )
            semi_min = self.pore_centre.distance( stoma.tracked_points[ 'pore-width' ][0] )

            pore_bdy_pts = [ _pt for _pt in pts if ( _pt.x / semi_maj )**2 + ( _pt.y / semi_min )**2 < 1 + padding ]

            # sort by angle - atan2 result is in (-pi,pi]
            pore_bdy_pts.sort( key=lambda _p: math.atan2( _p.y, _p.x ) )

            stoma.tracked_points[ 'pore-boundary' ] = pore_bdy_pts

        def create_calculators():
            """ set up the calculators for width, length etc. """

            # when only one guard cell is present we need to double the stoma and pore width because
            # it is measured from the pore centre
            double_if_1gc_fn = None if self.general_store.num_gcs == 2 else lambda x : 2 * x

            stoma_width_calc = xcalc.DistanceCalculator( prefix='stoma-width',
                                                         node_pair=stoma.tracked_points[ 'stoma-width' ],
                                                         lambda_fn=double_if_1gc_fn )

            stoma.calculators.append( stoma_width_calc )

            pore_width_calc = xcalc.DirectionalDistanceCalculator( prefix='pore-width',
                                                                   node_pair=stoma.tracked_points[ 'pore-width' ],
                                                                   direction=point.y_axis,
                                                                   lambda_fn=double_if_1gc_fn )

            stoma.calculators.append( pore_width_calc )

            pore_length_calc = xcalc.DistanceCalculator( prefix='pore-length',
                                                         node_pair=stoma.tracked_points[ 'pore-length' ] )

            stoma_length_calc = xcalc.DistanceCalculator( prefix='stoma-length',
                                                          node_pair=stoma.tracked_points[ 'stoma-length' ] )

            stoma.calculators.append( pore_length_calc )
            stoma.calculators.append( stoma_length_calc )

            # create calculator for the pore area
            pore_area_calc = xcalc.AreaCalculator2D( 'pore-area',
                                                     boundary_pts=stoma.tracked_points[ 'pore-boundary' ],
                                                     lambda_fn=double_if_1gc_fn )

            stoma.calculators.append( pore_area_calc )

            """ Store the height and width calculators """

            for gc_index in range( self.general_store.num_gcs ):
                for hw in ( 'height', 'width' ):
                    gc_key = 'gc{}-{}'.format( gc_index, hw )
                    pts = stoma.tracked_points[ gc_key ]

                    stoma.calculators.append( xcalc.DistanceCalculator( prefix=gc_key,
                                                                        node_pair=pts ) )

        print( '--> Searching for points of interest...' )

        padding = 0.001

        stoma = self.general_store.model_part_store.stoma

        # Store the tracked points
        #
        store_width_points()
        store_length_points()
        store_height_points()
        store_pore_boundary_nodes()

        self.general_store.model_part_store.set_tracked_point_list()

        create_calculators()

        if self.verbose:
            for mp in self.general_store.model_part_store.model_parts:
                mp.print_detail()

    def process_state_data( self, state_section ):
        """ Process the different parts of the state data for one state

        Example data:

        PLT_STATE                0x02000000  0x00101824
        PLT_STATE_HEADER         0x02010000  0x0000000c
        PLT_STATE_HDR_TIME       0x02010002  0x00000004  0x00000000
        PLT_STATE_DATA           0x02020000  0x00101808
        PLT_NODE_DATA            0x02020300  0x0001de98
        PLT_STATE_VARIABLE       0x02020001  0x0001de90
        PLT_STATE_VAR_ID         0x02020002  0x00000004  0x00000001
        PLT_STATE_VAR_DATA       0x02020003  0x0001de7c  0x00000000  0x0001de74  0x00000000 ...

        :param state_section: State object
        """

        def get_state_data_byte_offsets( byte_data ):
            """ Calculate the offsets of the bytes for the state data sections.
            :param byte_data:
            :return: dict
            """

            the_byte_offsets = dict( )

            for _id in hdr.FEBIO_STATE_DATA_IDS:
                _offset = xplt.get_offset_for_id( byte_data, _id, 0, force_dword=True )
                the_byte_offsets[ _id ] = None if _offset is None else _offset * hdr.SZ_DWORD

            return the_byte_offsets

        state_data = state_section.get_datum( )

        offset = 0

        state_header_section, offset = self.read_section( state_data, 'PLT_STATE_HEADER', offset )

        state_header_time_section, offset = self.read_section( state_data,
                                                               'PLT_STATE_HDR_TIME',
                                                               offset,
                                                               include_data=True )

        state_data_section, offset = self.read_section( state_data, 'PLT_STATE_DATA', offset )

        if self.verbose:
            print( state_section )
            print( state_header_section )
            print( state_header_time_section )
            print( state_data_section )

        # store the State object (time and attributes)
        state_obj = xcls.State( state_header_time_section.get_datum( ) )
        self.general_store.states.append( state_obj )

        print( '--> Processing state: {}. '.format( state_obj ), end='' )
        start_time = time.time()

        byte_offsets = get_state_data_byte_offsets( state_data )

        for state_data_id in hdr.FEBIO_STATE_DATA_IDS:
            offset = byte_offsets[ state_data_id ]

            if offset is None:
                continue

            # 'tag' will be PLT_xxx_DATA (xxx is NODE, ELEMENT etc.)
            tag = xplt.get_tag_for_id( state_data_id )

            state_data_section, offset = self.read_section( state_data, tag, offset, True )

            if self.verbose:
                print( state_data_section )

            self.process_state_section_data( state_obj, state_data_section )

        print( 'Time taken: {:5.3f}'.format( time.time() - start_time ) )

    @staticmethod
    def get_pt_from_dict( nid, nid_pt_dict, default=None ):
        """
        Get a point from a { node_id : point } dictionary for a given node id (aka nid)

        :param nid: required node_id
        :param nid_pt_dict: dict
        :param default: value to return if not found
        :return: the point - can be None
        """

        if nid is None:
            return default

        if nid in nid_pt_dict:
            return nid_pt_dict[ nid ]

        return default

    def process_state_section_data( self, state, state_data_section ):
        """ Process one of the state sections storing relevant information.
        The data consists of one or more sets of...

        PLT_STATE_VARIABLE       0x02020001  0x0001de90
        PLT_STATE_VAR_ID         0x02020002  0x00000004  0x00000001
        PLT_STATE_VAR_DATA       0x02020003  0x0001de7c  0x00000000  0x0001de74  0x00000000 ...

        N.B. First two data DWORDS are descriptive:         ^^^^        ^^^^
             The first one refers to a 'set' and the second is the size in bytes of remaining data

        :param state:
        :param state_data_section:
        """

        # Check the dictionary (belts and braces)
        #
        dict_list = self.general_store.model_dictionary.get( state_data_section.id, None )
        if dict_list is None:
            print( '** Found data for {} but no dictionary entry'.format( state_data_section.get_tag( ) ) )
            print( '** Dictionary is {}'.format( self.general_store.model_dictionary ) )
            return

        data = state_data_section.get_datum( )

        in_node_section = ( state_data_section.id == hdr.FEBIO_TAGS[ 'PLT_NODE_DATA' ] )

        st = self.general_store.model_part_store.stoma

        offset = 0
        while True:
            plt_state_var, offset = self.read_section( data, 'PLT_STATE_VARIABLE', offset )

            if plt_state_var is None:
                break

            plt_state_var_id, offset = self.read_section( data, 'PLT_STATE_VAR_ID', offset, True )
            plt_state_var_data, offset = self.read_section( data, 'PLT_STATE_VAR_DATA', offset, True )

            # get the dictionary item here so that we know how to process what's coming next
            # (-1 because lists are 0-based)
            dict_item = dict_list[ plt_state_var_id.get_datum_int( ) - 1 ]

            if self.verbose:
                print( plt_state_var )
                print( plt_state_var_id )
                print( plt_state_var_data )

                # first DWORD is a region id (0 indicates everything)
                # region_id = plt_state_var_data.get_datum( )[ 0 ]
                region_id = plt_state_var_data.get_datum_ndarray( )[ 0 ]
                print( 'Region id = {}'.format( region_id ) )

                print( 'Dictionary item: ', dict_item )

            # Could move this to post processing - will need to store displacements in the state object
            # although the advantage of doing it here is it does not need to be stored

            if in_node_section and dict_item.get_item_name( ) == 'displacement':
                raw_values = plt_state_var_data.get_datum_ndarray( )[ 2: ]

                xyz = np.frombuffer( raw_values, dtype=np.float32, offset=0 )
                displacements = self.xyz_list_to_dict( xyz )

                # Store the displacements for each tracked point
                #
                for tracked_pt in st.tracked_points[ 'tracked-points' ]:
                    displaced_pt = displacements[ tracked_pt.id ]
                    state.add_displacement( tracked_pt.id, displaced_pt )

                # Store all the displacements if we are computing the mesh quality
                if self.metrics.is_mesh_calculation_on:
                    for displaced_pt in displacements.values():
                        state.add_displacement( displaced_pt.id, displaced_pt )

                # Run the calculators...
                #
                displaced_tracked_pts = { pt.id: pt + displacements[ pt.id ]
                                          for pt in st.tracked_points[ 'tracked-points' ] }

                for calculator in st.calculators:
                    result = calculator.calculate( displaced_tracked_pts )

                    for key, value in result.items():
                        state.set_attribute( key, value )

                # TODO move into calculators (firstly, move the node_set into tracked points)...

                # Calculate the surface area and enclosed volume of each surface
                #
                for surface_idx, surface_data in enumerate( st.closed_surfaces.values() ):
                    # this happens when a dummy surface is added to the stoma due to the absence of the tip wall
                    if len( surface_data ) == 0:
                        continue

                    facets = surface_data.get( 'facets' )
                    node_set = surface_data.get( 'node_set' )

                    # create dict mapping node id to the displaced point
                    displaced_pts = { node_id: self.general_store.pts_t0.get( node_id ) + displacements[ node_id ]
                                      for node_id in node_set }

                    volume, area = geom.calculate_volume_and_area( displaced_pts, facets )

                    if self.verbose:
                        print( 'State: Area = {}, Volume = {}'.format( area, volume ) )

                    state.set_attribute( 'gc{}-surface-area'.format( surface_idx ), area )
                    state.set_attribute( 'gc{}-volume'.format( surface_idx ), volume )

    def read_file_to_string( self ):
        """ Open the XPLT file and read the whole thing (will close automatically) """

        print( '--> Reading file...' )

        with open( self.xplt_filename, 'rb' ) as fin:
            byte_data = np.fromfile( fin, dtype=np.uint8 )

        return byte_data

    @staticmethod
    def check_control_dword( data ):
        """ Make sure the first DWORD matches the FEBio specification
        :param data:
        """

        file_header = xplt.convert_4bytes_to_dword( data )

        if xplt.is_header_dword( file_header ):
            print( '--> The header DWORD is {:#010x}, which matches the expected value'.format( file_header ) )
        else:
            raise XpltReaderException(
                'The header DWORD is {:#010x}, which does not match the expected value'.format( file_header ) )

    def read_section( self, xplt_data, xplt_tag, offset, include_data=False ):
        """ Read a section from the raw data and update the offset.
        The first two DWORDS are the TAG and the SIZE

        :param xplt_data: the raw data (byte array)
        :type xplt_data: np.ndarray
        :param xplt_tag: the name of the FEBio tag
        :type xplt_tag: str
        :param offset: byte offset at which to start
        :type offset: int
        :param include_data: flag to include the data (optional, default=False) but may be ignored if
                the xplt_tag is associated with data (via FEBIO_TAGS_WITH_DATA)
        :type include_data: bool
        :return: a tuple containing the Section object and the new offset
        :rtype: ( xcls.Section, int )
        """

        if len( xplt_data ) <= offset:
            return None, None

        # Read the xplt id
        xplt_id = xplt.convert_4bytes_to_dword( xplt_data, offset )

        if xplt_id != hdr.FEBIO_TAGS[ xplt_tag ]:
            msg = 'Expected to find tag ({}) but found...'.format( xplt_tag )
            print( msg )
            self.output_raw_values( xplt_data, offset, 100 )

            raise XpltReaderException( msg )

        # move past the TAG
        offset += hdr.SZ_DWORD

        # the size is an array of 4-bytes so convert to a DWORD to get the size (which will be in bytes)
        sz_bytes = xplt.convert_4bytes_to_dword( xplt_data, offset )

        # move past the size
        offset += hdr.SZ_DWORD

        section_data = None
        if include_data or xplt_id in hdr.FEBIO_TAGS_WITH_DATA:
            section_data, offset = xplt.get_section_data( xplt_data, xplt_id, sz_bytes, offset )

        sect = xcls.Section( xplt_id, section_data )

        return sect, offset

    def do_post_processing( self ):
        """
        Perform post processing on the data obtained from the XPLT file
        :return:
        """

        print( OUTPUT_DELIMITER )
        print( '--> Performing post-processing...' )

        # calculate the aspect ratio (at the mid point) for each guard cell
        for state in self.general_store.states:
            for gc_key in ( 'gc0', 'gc1' ):
                key_width, key_height = [ '{}-{}'.format( gc_key, wh ) for wh in ( 'width', 'height' ) ]

                if state.has_attribute( key_width ) and state.has_attribute( key_height ):
                    width, height = state.get_attribute( key_width ), state.get_attribute( key_height )
                    state.set_attribute( '{}-ar'.format( gc_key ), height / width )

        if self.metrics.is_compare_vs_open_stoma_on:
            print( '--> Calculating differences vs. open stoma...' )

            for state in self.general_store.states:
                result = self.metrics.evaluate_metric( state )

                state.set_attribute( 'optimisation-keys', [ _[0] for _ in result ] )

                for kv_pair in result:
                    state.set_attribute( kv_pair[0], kv_pair[1] )

                # Store all of the diffs in one place (for the LM optimiser)
                state.set_attribute( 'all-diffs', [ _[1] for _ in result if 'metric' not in _[0] ] )

        # Shall we check the mesh quality?
        if self.metrics.is_mesh_calculation_on:
            node_ids = sorted( list( self.general_store.pts_t0.keys( ) ) )

            for state in [ self.general_store.states[ idx ] for idx in (0, -1) ]:
                print( '--> Calculating mesh metrics using VTK for state {} (t={})...'.format(
                    state.ordinal, state.time ) )

                disps = state.displacements

                pts = [ self.general_store.pts_t0[ nid ] + disps[nid ] for nid in node_ids ]

                m_c = mesh_ch.VTKMeshChecker.create_from_points_and_elements( pts,
                                                                              self.general_store.elements )

                mesh_quality = m_c.calculate_mesh_quality( )
                for mesh_metric in mesh_quality:
                    print( '--> {}'.format( mesh_metric ) )

                # write mesh to a vtk file
                fn_split = os.path.splitext( self.xplt_filename )
                m_c.output_vtk_file( file_name='{}.state={}.vtk'.format( fn_split[0], state.ordinal ) )

                state.set_attribute( 'mesh-quality', mesh_quality )

    def process_file( self, results_filename, output_raw=False ):
        """
        The XPLT file is a hierarchical structure of DWORDs (i.e. uint32).
        The PLT_ROOT section is first with header, material and geometry sections (to name a few)
        followed by a PLT_STATE section for each time point.

        :param results_filename:
        :param output_raw:
        """

        print( OUTPUT_DELIMITER )

        xplt_data = self.read_file_to_string( )

        self.check_control_dword( xplt_data )

        if output_raw:
            self.output_raw_values( xplt_data )
            return

        # Process the root section
        #
        print( '--> Processing root section...' )
        offset = 4

        root_section, offset = self.read_section( xplt_data, 'PLT_ROOT', offset, True )

        self.process_root_section( root_section )

        # Process the state sections...
        #
        print( '--> Processing state sections...' )

        # the size may have been altered (due to the DWORD realignment) so find PLT_STATE
        offset = xplt.get_offset_for_tag( xplt_data, 'PLT_STATE', offset )

        while True:
            state_section, offset = self.read_section( xplt_data, 'PLT_STATE', offset, True )

            if state_section is None:
                break

            self.process_state_data( state_section )

            # # Useful for debugging:
            # if len( self.general_store.states ) > 3:
            #     break

        print( '--> Processed {} state sections'.format( len( self.general_store.states ) ) )

        self.do_post_processing( )

        # save data
        self.write_results_to_file( results_filename )

        self.general_store.print_detail( )


    def write_results_to_file( self, results_filename ):
        """ Write the results to a file with a specific extension """

        # Is there any state data ?
        if len( self.general_store.states ) == 0:
            print( '--> No state data' )
            return

        rfw = ResultsFileWriter( self.general_store, self.metrics )

        rfw.write_to_file( results_filename )


def _format_results_filename( xplt_filename, results_filename ):

    if results_filename is None:
        # Create file for the processed results (get rid of the extension)
        f_name = os.path.splitext( xplt_filename )[ 0 ]

        results_filename = '{}.stats.txt'.format( f_name )

    if not results_filename.endswith( '.stats.txt' ):
        results_filename += '.stats.txt'

    return results_filename


def process_xplt_file( xplt_filename,
                       results_filename=None,
                       verbose=False,
                       output_raw_data=False,
                       metrics=None ):

    """ Extract useful information from an XPLT file.

    :param xplt_filename: the name of the XPLT file
    :param results_filename: name of the results file (will be set if omitted)
    :param verbose: verbosity flag
    :param output_raw_data: useful for debugging (produces a lot of output)
    :param metrics: the metrics to calculate :py:class:`xcalc.XpltReaderMetrics`
    :return: the final state from the simulation (or None if the simulation failed)
    :rtype: :py:class:`xcls.State`
    """

    if metrics is None:
        metrics = xcalc.XpltReaderMetrics()

    xplt_file_reader = XpltFileReader( xplt_filename, verbose, metrics )

    final_state = None

    if os.path.exists( xplt_filename ):
        results_filename = _format_results_filename( xplt_filename=xplt_filename, results_filename=results_filename )

        start_time = time.time( )
        xplt_file_reader.process_file( output_raw=output_raw_data, results_filename=results_filename )
        end_time = time.time( )

        print( '' )
        print( '--> Time to process file : {:5.3f} seconds'.format( end_time - start_time ) )
        print( '' )

        final_state = xplt_file_reader.general_store.states[ -1 ]
    else:
        print( '--> File does not exist' )

    return final_state


if __name__ == '__main__':
    raise NotImplementedError( 'Access to this functionality is via the top-level module' )
