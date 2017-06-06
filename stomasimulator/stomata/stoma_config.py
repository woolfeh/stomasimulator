""" Dimensions of a stoma """

from __future__ import print_function

from math import log, sqrt

from sortedcontainers import SortedDict

import stomasimulator.geom.ellipse as el
import stomasimulator.febio.feb.feb_material as fm

class PressureDefinition(object):
    def __init__(self, yaml_def):
        self.stoma_max = yaml_def[ 'stoma' ][ 'max' ]
        self.epidermis_max = yaml_def[ 'epidermis' ][ 'max' ] if 'epidermis' in yaml_def else 0.0


def generate_stoma_model( yaml_file, model_label, material_label, whole_stoma=False, optimise_vs_open=True ):
    """
    Reformat the dict obtained from the YAML-based configuration file so that it's suitable for the
      StomateConfig class

    :param yaml_file: representation of the contents of the stoma configuration file
    :type yaml_file: config.YAMLFileContents
    :param model_label: label of the model - it's the key under 'model -> dimensions' in cfg_dict
    :param material_label: label for the material model, e.g. tiVW
    :param optimise_vs_open: use the 'open' dimensions as the reference dimensions in the comparison helper
    :return: a StomateConfig object
    :rtype: StomaConfig
    """

    # TODO: replace this factory method with direct construction of StomateConfig object

    def get_model_sequence():
        """ Get a list of the models in the hierarchy with the root model first """
        sequence = [ model_label, ]

        while True:
            model_definition = yaml_file.models[ sequence[-1] ]
            if 'ref-model' in model_definition:
                sequence.append( model_definition[ 'ref-model' ] )
            else:
                break

        sequence.reverse()

        return sequence

    def get_material_model():
        # Create appropriate material model (default parameters will be set)
        #
        if material_label is None:
            return None

        mm = fm.lookup_material_model( material_model_label=material_label )

        if material_label in material_def:
            material_params = material_def[ material_label ]

            # Use the values to override the defaults

            if mm.is_isotropic:
                mm.set_parameters( c1=material_params[ 'c1' ],
                                   c2=material_params[ 'c2' ] )
            else:
                mm.set_parameters( c1=material_params[ 'c1' ],
                                   c2=material_params[ 'c2' ],
                                   c5=material_params[ 'c5' ],
                                   lm=material_params[ 'lm' ] )

        return mm


    open_geom = dict()
    closed_geom = dict()
    material_def = dict()
    pressure_def = dict( yaml_file.pressure )

    for model_label in get_model_sequence():
        model_def = yaml_file.models[ model_label ]

        # process the observation
        if 'observation' in model_def:
            observation_label = model_def[ 'observation' ]
            observation_def = yaml_file.observations[ observation_label ]

            open_geom.update( observation_def[ 'open' ] )
            closed_geom.update( observation_def[ 'closed' ] )

        if 'geometry' in model_def:
            closed_geom.update( model_def[ 'geometry' ] )

        if 'material' in model_def:
            material_def.update( model_def[ 'material' ] )

        if 'pressure' in model_def:
            pressure_def.update( model_def[ 'pressure' ] )

    wt = WallThicknesses( closed_geom )

    mm = get_material_model()

    ch = ComparisonHelper( reference_dimensions=open_geom if optimise_vs_open else closed_geom,
                           optimisation_def=yaml_file.optimisation,
                           open_pressure=yaml_file.pressure[ 'stoma' ][ 'max' ] )

    sc = StomaConfig( label=model_label,
                      whole_stoma=whole_stoma,
                      closed_geom=StomaGeometry.from_dict( closed_geom ),
                      open_geom=StomaGeometry.from_dict( open_geom, ignore_missing=True ),
                      material_model=mm,
                      wall_thickness=wt,
                      pressure=PressureDefinition( pressure_def ),
                      comparison_helper=ch )

    return sc


class ComparisonHelper(object):
    """ Helps comparisons vs. the observed open stoma """

    def __init__(self, reference_dimensions, optimisation_def, open_pressure ):
        """
        :param reference_dimensions: The observed dimensions
        :type reference_dimensions: dict
        :param optimisation_def: The keys and any aliases for the optimisation. The aliases provide mapping(s) for any
        names that differ
        :type optimisation_def: dict
        :param open_pressure: The guard cell turgor pressure when the stoma is fully open
        :type open_pressure: float
        """

        self.reference_dimensions = SortedDict( reference_dimensions )
        self.optimisation_keys = sorted( optimisation_def[ 'keys' ] )
        self.key_aliases = optimisation_def[ 'aliases' ]
        self.open_pressure = open_pressure

        # make sure the keys on which the optimisation will be performed are in the observed set
        if len( set( self.optimisation_keys ) - set( reference_dimensions.keys() ) ) > 0:
            raise KeyError( "One or more of the optimisation keys {} is/are missing from the reference dimensions {}".
                            format( set( self.optimisation_keys ), set( self.reference_dimensions.keys() ) ) )

        # make sure any key aliases for the optimisation are in the observed set
        if len( set( self.key_aliases.keys() ) - set( reference_dimensions.keys() ) ) > 0:
            raise KeyError( "One or more of the key aliases {} is/are missing from the reference dimensions {}".
                            format( set( self.key_aliases.keys() ), set( self.reference_dimensions.keys() ) ) )

    @property
    def simulation_keys(self):
        """
        Get the keys produced by the simulation - names are set in the configuration file
        :return: list
        """
        return [ self.key_aliases[ key ] if key in self.key_aliases else key for key in self.optimisation_keys ]

    def perform_comparison(self, state_pressure, state_data):
        """
        Calculate the metric and percent difference to each measurement

        :param state_pressure: pressure at the simulation state
        :type state_pressure: float
        :param state_data: the data extracted for the state
        :type state_data: dict

        :return: Each tuple comprises 2 items: a label and then its numerical value
        :rtype: list of tuple
        """

        ref_dimensions = [ self.reference_dimensions[ key ] for key in self.optimisation_keys ]
        sim_dimensions = [ state_data[ key ] for key in self.simulation_keys ]

        # the difference will be -1 if the simulation value is zero - happens to pore-width when the stoma closes
        pc_diffs = [ sim / ref - 1.0 for ref, sim in zip( ref_dimensions, sim_dimensions ) ]

        # metric is 'sum [ ln( predict_i / actual_i )**2 ]
        metric_raw = sum( log( max( abs( 1 + pc_diff ), 1e-5 ) ) ** 2 for pc_diff in pc_diffs )
        metric_raw = sqrt( metric_raw )

        dp = self.open_pressure - state_pressure

        # weight the metric to the end time/pressure
        metric_weighted = metric_raw + dp

        result_keys = [ 'metric', 'metric_raw', 'dp' ] + [ 'pc-{}'.format( k ) for k in self.simulation_keys ]
        result_vals = [ metric_weighted, metric_raw, dp ] + pc_diffs

        return zip( result_keys, result_vals )


class StomaGeometry(object):
    def __init__(self, pore_width, pore_length, gc_width, stoma_length, ignore_missing=False, **kwargs):
        self.pore_length = pore_length
        self.pore_width = pore_width
        self.gc_width = gc_width
        self.stoma_length = stoma_length

        # store the remaining values
        for key, value in kwargs.items():
            key_ = key.replace( '-', '_')
            self.__setattr__( key_, value )

        # Sanity check the parameters...
        #
        if not ignore_missing and self.pore_length and self.stoma_length:
            if self.pore_length > self.stoma_length:
                exit( "The pore length ({}) is greater than the stoma length ({}). Exiting...".
                      format( self.pore_length, self.stoma_length ) )

        if not ignore_missing and self.pore_width:
            if self.pore_width < 0.0:
                exit( 'The pore width ({}) must be positive. Exiting...'.format( self.pore_width ) )

            if self.stoma_width:
                if self.pore_width > self.stoma_width:
                    exit( "The pore width ({}) is greater than the stoma width ({}). Exiting...".
                          format( self.pore_width, self.stoma_width ) )

    @property
    def stoma_width(self):
        return self.pore_width + 2 * self.gc_width

    @staticmethod
    def from_dict(d, ignore_missing=False):
        """ Most of the time these are available but check anyway """

        pl = d[ 'pore-length' ]  if 'pore-length' in d else None
        pw = d[ 'pore-width' ]   if 'pore-width' in d else None
        sl = d[ 'stoma-length' ] if 'stoma-length' in d else None
        gcw = d[ 'gc-width' ]    if 'gc-width' in d else None

        d_ = dict( d )
        for key in ( 'pore-length', 'pore-width', 'stoma-length', 'gc-width' ):
            if key in d_:
                d_.pop( key )

        return StomaGeometry( pore_width=pw,
                              pore_length=pl,
                              gc_width=gcw,
                              stoma_length=sl,
                              ignore_missing=ignore_missing,
                              **d_ )


class WallThicknesses( object ):
    """ Store the wall thicknesses """

    # pylint: disable=too-few-public-methods

    def __init__( self, d ):
        """ Initialise from dictionary
        :param d:
        """
        self.dorsal = d[ 'cwt-dorsal' ]
        self.ventral = d[ 'cwt-ventral' ]
        self.periclinal = d[ 'cwt-periclinal' ]
        self.polar = d[ 'cwt-tip' ]
        self.epidermal = d[ 'cwt-epidermal' ] if 'cwt-epidermal' in d else None
        self.is_thicken_vp_on = 'cwt-thicken-vp' in d

    def print_detail( self ):
        """ Output information about the cell wall thicknesses """

        print( ' Cell wall thicknesses:' )
        print( '             Dorsal wall : {}'.format( self.dorsal ) )
        print( '         Periclinal wall : {}'.format( self.periclinal ) )
        print( '            Ventral wall : {}'.format( self.ventral ) )
        print( '        Thicken v-p wall : {}'.format( self.is_thicken_vp_on ) )
        print( '                Tip wall : {}'.format( self.polar ) )
        print( '          Epidermal wall : {}'.format( 'Not present' if self.epidermal is None
                                                       else self.epidermal ) )

        print( ' ' )


class StomaConfig( object ):
    """ The lengths, widths and cell wall thicknesses of a stoma """

    def __init__(self, label, whole_stoma, closed_geom, open_geom,
                 material_model, wall_thickness, pressure, comparison_helper=None ):
        """
        Initialise the stoma model definition

        :param label: str - label of the model
        :param whole_stoma: bool - 1 or 2 guard cells
        :param closed_geom: StomaGeometry - geometry of the closed stoma
        :param open_geom: StomaGeometry - geometry of the open stoma
        :param material: MaterialModel - definition of the material
        :param wall_thickness: WallThicknesses - the wall thicknesses
        :param comparison_helper: ComparisonHelper - performs calculations pertinent to the optimisation
        """
        self.label = label
        self.whole_stoma = whole_stoma
        self.closed = closed_geom
        self.open = open_geom
        self.material_model = material_model
        self.wall_thickness = wall_thickness
        self.pressure = pressure
        self.comparison_helper = comparison_helper

        # ellipse objects
        self.__ventral_ellipse = None
        self.__mid_ellipse = None
        self.__dorsal_ellipse = None

        self.__inner_ventral_curve = None
        self.__inner_dorsal_curve = None

        # Sanity check the parameters...
        #
        if self.mid_gc_width - self.wall_thickness.dorsal - self.wall_thickness.ventral < 0.0:
            exit( 'The GC width must be greater than the combined dorsal and ventral wall thicknesses. Exiting...' )


    def __str__( self ):
        return '[{cn}: stoma length/width={sl}/{sw}, pore length/width={pl}/{pw}]'.format(
            cn=self.__class__.__name__,
            sl=self.stoma_length, sw=self.stoma_width, pl=self.pore_length, pw=self.pore_width )

    @property
    def is_epidermal_wall_present( self ):
        return self.wall_thickness.epidermal is not None

    def print_detail( self ):
        """ Output information about the dimensions """

        title = "'{}' stoma configuration parameters".format( self.label )
        bdy = '=' * len( title )

        print( bdy )
        print( title )
        print( bdy )
        print( ' ' )
        print( ' Geometry parameters' )
        print( ' -------------------' )
        print( ' Stoma:   length / width : {} / {}'.format( self.stoma_length, self.stoma_width ) )
        print( ' Pore:    length / width : {} / {}'.format( self.pore_length, self.pore_width ) )
        print( ' GC: tip length / height : {} / {}'.format( self.tip_length, self.tip_height ) )
        print( ' GC:  mid width / height : {} / {}'.format( self.mid_gc_width, self.mid_gc_height ) )
        print( ' GC:        aspect ratio : {} (>1 => bulges out of leaf plane)'.format( self.aspect_ratio ) )
        print( '' )
        print( ' Number of guard cells   : {}'.format( self.num_gcs ) )
        print( '' )

        self.wall_thickness.print_detail()

    @property
    def pore_length(self):
        return self.closed.pore_length

    @property
    def pore_width(self):
        return self.closed.pore_width

    @property
    def gc_width(self):
        return self.closed.gc_width

    @property
    def stoma_length(self):
        return self.closed.stoma_length

    @property
    def stoma_width(self):
        return self.closed.stoma_width

    @property
    def aspect_ratio(self):
        return self.closed.aspect_ratio

    @property
    def num_gcs(self):
        return 2 if self.whole_stoma else 1

    @property
    def mid_gc_width( self ):
        """ Width of the GC at the mid point (in the equatorial plane) """
        return (self.stoma_width - self.pore_width) / 2

    @property
    def mid_gc_height( self ):
        """ Height at the GC mid point """
        return self.mid_gc_width * self.aspect_ratio

    @property
    def tip_length( self ):
        """ Length of the tip (where the GCs meet) in the equatorial plane """
        return (self.stoma_length - self.pore_length) / 2

    @property
    def tip_height( self ):
        """ Height of the tip (where the GCs meet) """
        return self.tip_length * self.aspect_ratio

    # for the ellipses (use a,b etc. for shorthand convenience)...

    @property
    def dorsal_ellipse( self ):
        """ Ellipse that approximates the dorsal wall """
        if self.__dorsal_ellipse is None:
            a_d = self.stoma_length / 2
            b_d = self.stoma_width / 2
            self.__dorsal_ellipse = el.Ellipse( a_d, b_d )

        return self.__dorsal_ellipse

    @property
    def inner_dorsal_curve( self ):
        """ Polygonal line that is offset from the dorsal wall by the wall thickness """

        if self.__inner_dorsal_curve is None:
            self.__inner_dorsal_curve = el.calculate_polyline_offset_from_ellipse( self.dorsal_ellipse,
                                                                                   -self.wall_thickness.dorsal )

        return self.__inner_dorsal_curve

    @property
    def mid_ellipse( self ):
        """ Ellipse that bisects the dorsal and ventral walls """
        if self.__mid_ellipse is None:
            a_mid = ( self.pore_length + self.tip_length  ) / 2
            b_mid = ( self.pore_width + self.mid_gc_width ) / 2

            self.__mid_ellipse = el.Ellipse( a_mid, b_mid )

        return self.__mid_ellipse

    @property
    def ventral_ellipse( self ):
        """ Ellipse that approximates the ventral wall """
        if self.__ventral_ellipse is None:
            a_v = self.pore_length / 2
            b_v = self.pore_width / 2

            self.__ventral_ellipse = el.Ellipse( a_v, b_v )

        return self.__ventral_ellipse

    @property
    def inner_ventral_curve( self ):
        """ Polygonal line that is offset from the ventral wall by the wall thickness """
        if self.__inner_ventral_curve is None:
            self.__inner_ventral_curve = el.calculate_polyline_offset_from_ellipse( self.ventral_ellipse,
                                                                                    self.wall_thickness.ventral )

        return self.__inner_ventral_curve

if __name__ == '__main__':
    raise RuntimeError( 'Access to this functionality is via the top-level module' )
