
import stomasimulator.febio.feb.feb_file_classes as ffc

class TestMethods(object):
    def test_get_timestep_control_1(self):
        tc = ffc.get_timestep_control( end_time=5.0, timestep=0.1 )

        assert tc.end_time == 5.0 and tc.dt == 0.1

    def test_get_timestep_control_1( self ):
        tc = ffc.get_timestep_control( end_time=1.0, timestep=0.2 )

        assert tc.end_time == 1.0 and tc.dt == 0.2

class TestLoadCurve(object):
    def test_lc_0(self):
        ffc.LoadCurve.reset_load_curve_id( )

        assert ffc.LoadCurve._load_curve_id == 1

    def test_lc_1(self):
        ffc.LoadCurve.reset_load_curve_id( )

        lc = ffc.LoadCurve()

        assert lc.load_curve_id == 1 and ffc.LoadCurve._load_curve_id == 2

    def test_lc_2(self):
        ffc.LoadCurve.reset_load_curve_id( )

        lc_1 = ffc.LoadCurve()
        lc_2 = ffc.LoadCurve()

        assert lc_1.load_curve_id == 1 and lc_2.load_curve_id == 2 and ffc.LoadCurve._load_curve_id == 3

    def test_lc_3( self ):
        ffc.LoadCurve.reset_load_curve_id( )

        lc_1 = ffc.LoadCurve()
        lc_2 = ffc.LoadCurve()

        ffc.LoadCurve.reset_load_curve_id( )

        assert lc_1.load_curve_id == 1 and lc_2.load_curve_id == 2 and ffc.LoadCurve._load_curve_id == 1

    def test_lc_4( self ):
        ffc.LoadCurve.reset_load_curve_id( )

        lc = ffc.LoadCurve( )
        lc.add_point( 0.0, 0.0 )

        assert len( lc.values ) == 1 and lc.values[0] == ( 0.0, 0.0 )

    def test_lc_5( self ):
        ffc.LoadCurve.reset_load_curve_id( )

        lc = ffc.LoadCurve( )
        lc.add_point( 0.0, 0.0 )
        lc.add_point( 1.0, 2.0 )

        assert len( lc.values ) == 2 and lc.values[0] == ( 0.0, 0.0 ) and lc.values[1] == ( 1.0, 2.0 )
