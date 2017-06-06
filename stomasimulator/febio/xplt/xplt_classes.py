""" Helper classes used when processing the XPLT file """

from __future__ import print_function

import collections

import numpy

import stomasimulator.geom.geom_utils as geom
import stomasimulator.stomata.stoma_config as sc

from . import xplt_utils as xut


class Section( object ):
    """ Represent a section of an XPLT file """

    def __init__( self, xplt_id, datum=None ):
        self.id = xplt_id
        self.datum = datum
        self.tag = xut.get_tag_for_id( self.id )

    def __str__( self ):
        nm = 'Unrecognised' if self.tag is None else self.tag

        fmt_str = '{:25s} {:#010x}'.format( nm, self.id )

        if self.datum is not None:
            fmt_str = '{} : datum {}'.format( fmt_str, self.datum )

        return fmt_str

    def has_tag( self ):
        """ Does the section have a FEBio tag """

        return self.tag is not None

    def get_tag( self, verbose=True ):
        """ Return the FEBio tag """

        if verbose and self.tag is None:
            print( 'Returning empty TAG for id : 0x%010x' % self.id )
        return self.tag

    def get_datum( self ):
        """
        :return: the datum (usually an ndarray)
        :rtype: list
        """
        return self.datum

    def get_datum_int( self ):
        """ Return datum converted to an int
        :return:
        :rtype: int
        """
        return int( self.datum )

    def get_datum_str( self, default='' ):
        """ Return datum converted to a string """
        return xut.convert_byte_array_to_str( self.datum, default=default )

    def get_datum_ndarray( self ):
        """
        :return:
        :rtype: numpy.ndarray
        """

        if isinstance( self.datum, numpy.ndarray ):
            return self.datum

        it = iter( self.datum )

        return numpy.ndarray( (len( it ),), buffer=self.datum )


class DictionaryItem( object ):
    """ Represent a dictionary item from an XPLT file """

    def __init__( self, item_type, item_fmt, item_name ):
        """ Initialise the dictionary item

        :param item_type:
        :type item_type: int
        :param item_fmt:
        :type item_fmt: int
        :param item_name:
        :type item_name: numpy.ndarray
        """
        self.item_type = item_type
        self.item_fmt = item_fmt
        self.item_name = xut.convert_byte_array_to_str( item_name, default='<Unknown>' )

    def __str__( self ):
        return 'DictItem: type={} ({}), format={} ({}), name=\'{}\''.format(
            self.get_item_type_str( ), self.item_type, self.get_item_format_str( ), self.item_fmt, self.item_name )

    def __repr__( self ):
        return '[' + str( self ) + ']'

    def get_item_type( self ):
        """ Return the item type
        :return:
        :rtype: int
        """
        return self.item_type

    def get_item_type_str( self ):
        """ Get the string representation of the item type
        :return:
        :rtype: str
        """
        if self.item_type == 0:
            return 'FLOAT'
        elif self.item_type == 1:
            return 'VEC3F'  # 3D vector of SP floats
        elif self.item_type == 2:
            return 'MAT3FS'  # symmetric 2nd order tensor of SP floats

        return '<Unknown type>'

    def get_item_format( self ):
        """ Get the item format
        :return:
        :rtype: int
        """
        return self.item_fmt

    def get_item_format_str( self ):
        """ Get a string representation of the format
        :return:
        :rtype: str
        """
        if self.item_type == 0:
            return 'NODE'  # one value for each node
        elif self.item_type == 1:
            return 'ITEM'  # one value for each item
        elif self.item_type == 2:
            return 'MULT'  # one value for each node for each item

        return '<Unknown type>'

    def get_item_name( self ):
        """ Get the name of the item """
        return self.item_name


class Material( object ):
    """ Represent a material """

    # pylint: disable=too-few-public-methods

    def __init__( self, mat_id, mat_name ):
        self.mat_id = mat_id
        self.mat_name = xut.convert_byte_array_to_str( mat_name, default='<unknown>' )

    def __str__( self ):
        return 'Material: id={}, name="{}" '.format( self.mat_id, self.mat_name )

    def __repr__( self ):
        return '[' + str( self ) + ']'


class State( object ):
    """ Represent a state in the XPLT file (time-point and data) """

    _ordinal = 1

    def __init__( self, time ):
        self.time = time
        self.ordinal = State._ordinal
        self.displacements = dict( )
        self.attributes = collections.OrderedDict( )

        State._ordinal += 1

    def __str__( self ):
        attr_str = "" if len( self.attributes ) == 0 else ', attributes={}'.format( str( self.attributes ) )
        return '[State {:4}, time={:7.6f}{}]'.format( self.ordinal, self.time, attr_str )

    def __repr__( self ):
        return str( self )

    def add_displacement( self, node_id, pt ):
        """ Add a displacement for a point
        :param node_id:
        :param pt:
        :return:
        """
        self.displacements[ node_id ] = pt

    def set_attribute( self, name, value ):
        """ Store area, volume, and anything else
        :param name:
        :param value:
        :return:
        """

        if isinstance( value, str ):
            # make sure the stored data does not contain spaces because
            # it could mess with column numbering in an output file
            self.attributes[ name ] = value.replace( ' ', '_' )
        else:
            self.attributes[ name ] = value

    def has_attribute( self, name ):
        """ Does the state have a specific attribute
        :param name:
        :return:
        """
        return name in self.attributes.keys( )

    def get_attribute( self, name ):
        """ Get an attribute by name
        :param name:
        :return:
        """
        return self.attributes[ name ]


if __name__ == '__main__':
    pass
