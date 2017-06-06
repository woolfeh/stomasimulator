
import stomasimulator.utils.utils as ut

class TestMethods( object ):
    ## Test normal vector computations

    # can't compute normal if points are collinear
    def test_module_exists( self ):
        assert ut.module_exists( 'sys' )

    def test_module_does_not_exist( self ):
        import sys

        bad_module = 'a' * 10
        while bad_module in sys.modules.keys():
            bad_module *= 10

        assert not ut.module_exists( bad_module )

    def test_get_datetime_now(self):
        import datetime as dt
        assert isinstance( ut.get_datetime_now(), dt.datetime )

    def test_get_elapsed_millis(self):
        dt1 = ut.get_datetime_now()

        # waste time
        r = range(1000000)
        s = sum(r)

        em = ut.get_elapsed_millis( dt1 )

        assert em > 0.0

    def test_get_diff_millis(self):
        dt1 = ut.get_datetime_now()

        # waste time
        r = range(1000000)
        s = sum(r)

        dt2 = ut.get_datetime_now()

        assert ut.get_diff_millis( dt1, dt2 ) > 0.0

    def test_convert_timedelta_to_millis(self):
        dt1 = ut.get_datetime_now()

        # waste time
        r = range(1000000)
        s = sum(r)

        dt2 = ut.get_datetime_now()

        td = dt2 - dt1

        assert ut.convert_timedelta_to_millis( td ) > 0.0
