""" Utility functions for handling the FEBio XPLT files """

from __future__ import print_function

import numpy as np

from . import xplt_header as hdr


def is_header_dword( a_dword ):
    """
    Does the given DWORD match FEBio's header DWORD?

    :param a_dword:
    :return: bool
    """
    return a_dword == hdr.PLT_FEBIO


def get_section_data( xplt_data, xplt_id, sz_bytes, offset ):
    """
    Read the data for a section

    :param xplt_data:
    :type xplt_data: np.ndarray
    :param xplt_id:
    :param sz_bytes:
    :param offset:
    :return: tuple of ( section's data, offset positioned after the data )
    """

    # These fields will need to be converted to DWORD arrays
    i_wanna_be_a_dword_array = (hdr.FEBIO_TAGS[ 'PLT_NODE_COORDS' ],
                                hdr.FEBIO_TAGS[ 'PLT_ELEMENT' ],
                                hdr.FEBIO_TAGS[ 'PLT_FACE' ],
                                hdr.FEBIO_TAGS[ 'PLT_STATE_VAR_DATA' ])

    # Check the length
    assert len( xplt_data ) + offset >= sz_bytes

    data_bytes = xplt_data[ offset: offset + sz_bytes ]

    if xplt_id in i_wanna_be_a_dword_array:
        # convert byte array to DWORD (aka 'u4') array
        assert sz_bytes % hdr.SZ_DWORD == 0
        section_data = np.frombuffer( data_bytes, 'u4' )
    elif sz_bytes == 4:
        if xplt_id == hdr.FEBIO_TAGS[ 'PLT_STATE_HDR_TIME' ]:
            section_data = convert_4bytes_to_float( data_bytes )
        else:
            section_data = convert_4bytes_to_dword( data_bytes )
    else:
        section_data = data_bytes

    return section_data, offset + sz_bytes


def get_offset_for_tag( data, tag_str, offset=0 ):
    """
    Get the offset for a given tag

    :param data:
    :type data: np.ndarray
    :param tag_str:
    :param offset:
    :return: return the offset
    """
    return get_offset_for_id( data, hdr.FEBIO_TAGS[ tag_str ], offset )


def get_tag_for_id( xplt_id ):
    """
    Get the tag string for a given tag id

    :param xplt_id: tag id
    :return: tag str or None if one is not found
    """
    import sys

    tag = None
    try:
        tag = hdr.TAGS[ xplt_id ]
    except KeyError:
        sys.exc_clear( )

    return tag


def get_offset_for_id( byte_data, id_, offset=0, verbose=False, force_dword=False ):
    """
    Search for the given id in the data

    :param byte_data:
    :param id_:
    :param offset:
    :param verbose:
    :param force_dword:
    :return:
    """
    exceptions = (hdr.FEBIO_TAGS[ 'PLT_DOM_ELEM_LIST' ], hdr.FEBIO_TAGS[ 'PLT_FACE_LIST' ])

    # These will be found and their offsets returned but their sub-elements will not
    # because the large number of them slows things down
    #
    # PLT_DOM_ELEM_LIST        0x01042200  0x000e3918
    # PLT_ELEMENT              0x01042201  0x00000014  0x00000001  ...
    # < loads of PLT_ELEMENTs >
    # ...
    #
    # PLT_FACE_LIST            0x01043200  0x00044390
    # PLT_FACE                 0x01043201  0x00000020  0x00000001  ...
    # < loads of PLT_FACEs >
    # ...

    def convert_dword_to_4bytes( a_dword ):
        """ Convert DWORD to a 4 byte array """
        return np.array( [ (a_dword & (0xff << pos)) >> pos for pos in (0, 8, 16, 24) ], 'u1' )

    if force_dword:
        found = False

        dword_data = np.frombuffer( byte_data, 'u4' )
        sz = len( dword_data )

        while offset < sz:
            dword = dword_data[ offset ]

            if dword == id_:
                found = True
                break

            # advance past the tag DWORD
            offset += 1

            if dword in hdr.TAGS:
                # Read the size and move past it
                dsz = dword_data[ offset ]
                offset += 1

                if dword in exceptions or dword in hdr.FEBIO_TAGS_WITH_DATA:
                    # Move past the data
                    offset += dsz / hdr.SZ_DWORD

        if verbose:
            print( 'Searching for 0x{:08x} {:25s}, offset = {}, found = {}'
                   .format( id_, hdr.TAGS[ id_ ], offset, found ) )

        return offset if found else None
    else:
        id_bytes = convert_dword_to_4bytes( id_ )
        offset = get_offset_for_bytes( byte_data, id_bytes, offset )
        return offset


def get_offset_for_bytes( byte_data, the_bytes, offset=0 ):
    """
    Find the offset for a series of bytes

    :param byte_data:
    :param the_bytes:
    :param offset:
    :return:
    """
    sz = len( byte_data )
    sz_tb = len( the_bytes )

    safety_offset = sz - sz_tb

    found = False

    while offset < safety_offset:

        while offset < safety_offset and the_bytes[ 0 ] != byte_data[ offset ]:
            offset += 1

        some_bytes = byte_data[ offset:offset + sz_tb ]

        if np.array_equal( the_bytes, some_bytes ):
            found = True
            break

        offset += 1

    return offset if found else None


def get_byte_array_for_tag( tag_str ):
    """
    Translate the tag to its corresponding integer value and then into its byte array representation
    :param tag_str:
    :return:
    """
    tag_id = hdr.FEBIO_TAGS[ tag_str ]

    # frombuffer needs an object that exposes the buffer interface
    buf = np.array( tag_id, dtype='u4' )

    return np.frombuffer( buf, 'u1' )


def convert_4bytes( byte_array, type_str, index=0 ):
    """
    Convert the 4 bytes from 'index' to the specified type

    :param byte_array:
    :param type_str: e.g. 'u1', 'u2', 'u4', 'f'
    :param index: (optional)
    :return:
    """
    assert len( byte_array ) >= 4 + index
    # ensure it's a numpy array so that it implements the buffer interface (reqd by 'frombuffer')
    buf = np.array( byte_array[ index: index + hdr.SZ_DWORD ] )
    return np.frombuffer( buf, type_str )[ 0 ]


def convert_4bytes_to_dword( byte_array, index=0 ):
    """
    Convert first 4 bytes of the array to a DWORD

    :param byte_array:
    :param index:
    :return:
    """
    return convert_4bytes( byte_array, 'u4', index )


def convert_4bytes_to_float( byte_array, index=0 ):
    """
    Convert first 4 bytes of the array to a float

    :param byte_array:
    :param index:
    :return:
    """
    return convert_4bytes( byte_array, 'f', index )


def convert_byte_array_to_str( data, default='' ):
    """
    Convert byte array (terminated by a null byte) to a string

    :param data:
    :param default:
    :return:
    """
    data_str = bytearray( data ).split( '\0' )

    if len( data_str ) == 0:
        return default

    # make sure it's a string
    the_str = str( data_str[ 0 ] )

    # check it's readable
    if len( the_str ) == 0 or (len( the_str ) > 0 and not the_str[ 0 ].isalnum( )):
        return default

    return the_str.strip( )
