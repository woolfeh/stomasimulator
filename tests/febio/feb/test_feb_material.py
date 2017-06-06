
import stomasimulator.febio.feb.feb_material as fm

c1 = 2.0
c2 = 3.0
c5 = 10.0
lm = 1.1

class TestIsotropicMatrix(object):

    def test_1(self):
        im = fm.IsotropicMatrix( c1=c1, c2=c2 )

        assert im.c1 == c1 and im.c2 == c2

    def test_mr_1(self):
        mr = fm.MooneyRivlinMatrix( c1=c1, c2=c2 )

        assert mr.initial_shear_modulus() == 2.0 * ( c1 + c2 )

    def test_vw_1( self ):
        mr = fm.VerondaWestmannMatrix( c1=c1, c2=c2 )

        assert mr.initial_shear_modulus( ) == c1 * c2

class TestFibres( object ):
    def test_1(self):
        ff = fm.Fibres( c5=c5, lm=lm )

        assert ff.c5 == c5 and ff.lm == lm

    def test_2(self):
        ff = fm.Fibres()

        ff.set_parameter_values( c5=c5, lm=lm )

        assert ff.c5 == c5 and ff.lm == lm

class TestMaterialModels(object):
    def test_mr(self):
        mm = fm.MooneyRivlinMaterialModel( c1=c1, c2=c2 )

        assert mm.isotropic_matrix.c1 == c1 and mm.isotropic_matrix.c2 == c2 and \
               mm.febio_name == 'Mooney-Rivlin' and mm.label == 'MR'

    def test_vw( self ):
        mm = fm.VerondaWestmannMaterialModel( c1=c1, c2=c2 )

        assert mm.isotropic_matrix.c1 == c1 and mm.isotropic_matrix.c2 == c2 and \
               mm.febio_name == 'Veronda-Westmann' and mm.label == 'VW'

    def test_timr( self ):
        mm = fm.TransIsoMooneyRivlinMaterialModel( c1=c1, c2=c2, c5=c5, lm=lm )

        assert mm.isotropic_matrix.c1 == c1 and mm.isotropic_matrix.c2 == c2 and \
               mm.fibres.c5 == c5 and mm.fibres.lm == lm and \
               mm.febio_name == 'trans iso Mooney-Rivlin' and mm.label == 'tiMR'


    def test_tivw( self ):
        mm = fm.TransIsoVerondaWestmannMaterialModel( c1=c1, c2=c2, c5=c5, lm=lm )

        assert mm.isotropic_matrix.c1 == c1 and mm.isotropic_matrix.c2 == c2 and \
               mm.fibres.c5 == c5 and mm.fibres.lm == lm and \
               mm.febio_name == 'trans iso Veronda-Westmann' and mm.label == 'tiVW'

class TestLookupMaterialModel(object):
    def test_lookup_material_model(self):
        mr = fm.lookup_material_model( 'MR' )
        vw = fm.lookup_material_model( 'VW' )
        timr = fm.lookup_material_model( 'tiMR' )
        tivw = fm.lookup_material_model( 'tiVW' )

        assert isinstance( mr, fm.MooneyRivlinMaterialModel ) and \
            isinstance( vw, fm.VerondaWestmannMaterialModel ) and \
            isinstance( timr, fm.TransIsoMooneyRivlinMaterialModel ) and \
            isinstance( tivw, fm.TransIsoVerondaWestmannMaterialModel )

