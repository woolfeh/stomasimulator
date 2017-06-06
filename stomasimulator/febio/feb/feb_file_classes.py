""" Some useful methods for the FEBio simulation """

from __future__ import division

import collections as coll

# Time step control parameters
TimestepParameters = coll.namedtuple( 'TimestepParameters',
                                      'num_steps  dt  end_time  min_dt  max_dt' )


def get_timestep_control( end_time, timestep=0.1 ):
    """
    Calculate the time stepping parameters for the FEBio simulation
    :param end_time:
    :param timestep:
    :return:
    """

    num_steps = int( end_time / timestep )
    min_dt = timestep / 10
    max_dt = timestep

    return TimestepParameters( num_steps=num_steps, dt=timestep, end_time=end_time,
                               min_dt=min_dt, max_dt=max_dt )


class LoadCurve( object ):
    """ Representation of a LoadCurve """

    _load_curve_id = 1

    def __init__( self ):
        self.load_curve_id = LoadCurve._load_curve_id
        self.values = list()

        LoadCurve._load_curve_id += 1

    def add_point(self, pt_time, value):
        self.values.append( ( pt_time, value ) )

    @staticmethod
    def reset_load_curve_id( ):
        """ Reset the curve id so that it counts from 1 """
        LoadCurve._load_curve_id = 1


if __name__ == '__main__':
    raise NotImplementedError( 'Access to this functionality is via the top-level module' )
