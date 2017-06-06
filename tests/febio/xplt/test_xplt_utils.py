##########################################################################################

import pytest
import numpy as np

import stomasimulator.febio.xplt.xplt_utils as x

##########################################################################################

class TestMethods( object ):
    """ Test the methods in hw_xplt """

    arr_3b = np.array( [1, 2, 3], dtype=np.uint8 )
    arr_4b = np.array( [1, 2, 3, 4], dtype=np.uint8 )
    arr_5b = np.array( [1, 2, 3, 4, 5], dtype=np.uint8 )

    # convert --> DWORD

    def test_convert_4bytes_to_dword_1( self ):
        with pytest.raises( AssertionError ):
            x.convert_4bytes_to_dword( TestMethods.arr_3b )

    def test_convert_4bytes_to_dword_2( self ):
        dword = x.convert_4bytes_to_dword( TestMethods.arr_4b )

        assert dword == 0x04030201

    def test_convert_4bytes_to_dword_3( self ):
        with pytest.raises( AssertionError ):
            x.convert_4bytes_to_dword( TestMethods.arr_4b, index=1 )

    def test_convert_4bytes_to_dword_4( self ):
        dword = x.convert_4bytes_to_dword( TestMethods.arr_5b )

        assert dword == 0x04030201

    def test_convert_4bytes_to_dword_5( self ):
        dword = x.convert_4bytes_to_dword( TestMethods.arr_5b, index=1 )

        assert dword == 0x05040302

    def test_convert_4bytes_to_dword_6( self ):
        with pytest.raises( AssertionError ):
            x.convert_4bytes_to_dword( TestMethods.arr_5b, index=2 )

    # convert --> float

    def test_convert_4bytes_to_float_1( self ):
        with pytest.raises( AssertionError ):
            x.convert_4bytes_to_float( TestMethods.arr_3b )

    def test_convert_4bytes_to_float_2( self ):
        a_float = x.convert_4bytes_to_float( TestMethods.arr_4b )

        np.testing.assert_allclose( a_float, 1.5399896e-36 )

    def test_convert_4bytes_to_float_3( self ):
        with pytest.raises( AssertionError ):
            x.convert_4bytes_to_float( TestMethods.arr_4b, 1 )

    def test_convert_4bytes_to_float_4( self ):
        a_float = x.convert_4bytes_to_float( TestMethods.arr_5b )

        np.testing.assert_allclose( a_float, 1.5399896e-36 )

    def test_convert_4bytes_to_float_5( self ):
        a_float = x.convert_4bytes_to_float( TestMethods.arr_5b, index=1 )

        np.testing.assert_allclose( a_float, 6.2071626e-36 )

    def test_convert_4bytes_to_float_6( self ):
        with pytest.raises( AssertionError ):
            x.convert_4bytes_to_float( TestMethods.arr_5b, index=2 )

##########################################################################################
