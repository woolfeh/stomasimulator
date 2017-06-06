"""
IO functions for datasets
"""

import os
import os.path as osp


class DataSet( object ):
    """
    Handle a dataset
    """

    def __init__( self ):
        """
        :return:
        """
        self.comments = list( )
        self.data = list( )

    def is_empty( self ):
        """
        Is there any data or comments
        :return: True if there are no comments and no data
        """
        return len( self.comments ) + len( self.data ) == 0

    def add_data_line( self, data_line ):
        """
        Add a data line
        :param data_line:
        :return:
        """
        self.data.append( data_line.strip( ) )

    def add_comment_line( self, comment_line ):
        """
        Add a comment line
        :param comment_line:
        :return:
        """
        self.comments.append( comment_line.strip( ) )

    def __repr__( self ):
        return '{} comment lines and {} data lines'.format( len( self.comments ), len( self.data ) )


class DataFile( object ):
    """
    Represent a datafile
    """

    def __init__( self, fname, datasets ):
        self._file_name = fname
        self._datasets = datasets

    @property
    def file_name( self ):
        """
        :return: file name
        """
        return self._file_name

    @property
    def datasets( self ):
        """
        list of dataset
        :return:
        """
        return self._datasets

    def __repr__( self ):
        return '{} - {} data sets'.format( self.file_name, len( self.datasets ) )

    def __iter__( self ):
        return self.datasets


class _DataFileParser( object ):
    """
    Data file parser
    """

    def __init__( self, fname ):
        self._full_file_name = fname

    @property
    def full_file_name( self ):
        """
        Full file name
        :return: string
        """
        return self._full_file_name

    def __iter__( self ):
        a_record = None
        l0 = ''

        with open( self.full_file_name, 'r' ) as fh:
            for line in fh:
                l1 = line.strip( )

                if len( l0 ) + len( l1 ) == 0 and a_record is not None:
                    yield a_record
                    a_record = None
                else:
                    if a_record is None:
                        a_record = DataSet( )

                    if len( l1 ) > 0:
                        if l1.startswith( '#' ):
                            a_record.add_comment_line( l1 )
                        else:
                            a_record.add_data_line( l1 )

                l0 = l1

        if a_record is not None:
            yield a_record


def read_files( dir_name ):
    """
    Read a set of files that end in 'stats.txt'

    :param dir_name:
    :return: a list of files
    :rtype: list
    """

    file_names = sorted( [ f for f in os.listdir( dir_name )
                           if osp.isfile( osp.join( dir_name, f ) ) and f.endswith( 'stats.txt' ) ] )

    data_files = list( )

    for fn in file_names:
        data_file = read_file( dir_name, fn )
        data_files.append( data_file )

    return data_files


def read_file( dir_name, file_name ):
    """
    Read a file into a DataFile object

    :param dir_name:
    :param file_name:
    :return: a datafile
    :rtype: DataFile
    """

    full_fname = osp.join( dir_name, file_name )

    data_sets = list( )
    for ds in _DataFileParser( full_fname ):
        data_sets.append( ds )

    df = DataFile( file_name, data_sets )

    return df


def get_time_course_for_columns( time_course_data_set, col_indices=None ):
    """
    Pick out a set of columns

    :param time_course_data_set:
    :param col_indices:
    :return:
    """

    import numpy as np

    if col_indices is None or None in col_indices:
        return None

    data = np.ndarray( shape=(len( time_course_data_set.data ), len( col_indices )), dtype=float )
    data.fill( None )

    for line_idx, data_line in enumerate( time_course_data_set.data ):
        split_data = data_line.split( )

        for idx, col_idx in enumerate( col_indices ):
            if 0 <= col_idx <= len( split_data ):
                data[ line_idx, idx ] = float( split_data[ col_idx ] )

    return data
