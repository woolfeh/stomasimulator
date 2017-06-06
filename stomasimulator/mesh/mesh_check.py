""" Call VTK to obtain mesh checking metrics """

from __future__ import print_function

import math

import stomasimulator.geom.point as gp

try:
    import vtk
except ImportError:
    # silence the error - vtk unavailable on the cluster
    # use an empty file to set the import reference
    from . import vtk_empty as vtk


class VTKMeshChecker(object):
    """ Handles mesh conversion and usage of the metrics defined in VTK """

    def __init__(self, vtk_mesh):
        """
        Create a mesh checker that uses the functionality in VTK.

        :param vtk_mesh: A VTK mesh

        :return:
        """

        self.vtk_mesh = vtk_mesh

    @classmethod
    def create_from_points_and_elements( cls, points, elements ):
        """
        Converts a list of points and a list of elements to a mesh checker

        :param points: List of Point objects
        :param elements: List of node ids for each element (correspond to the node id of the Point object)

        :return: Mesh checker
        :rtype: VTKMeshChecker
        """
        # create VTK points
        vtk_points = vtk.vtkPoints()

        for pt in points:
            vtk_points.InsertNextPoint( pt.xyz )

        # create the VTK mesh and add the points
        vtk_mesh = vtk.vtkUnstructuredGrid()
        vtk_mesh.SetPoints( vtk_points )

        # create the VTK elements
        for ele in elements:
            # NOTE: only hexahedra have metrics in VTK
            vtk_ele = vtk.vtkHexahedron() if len( ele ) == 8 else ( vtk.vtkWedge() if len( ele ) == 6 else None )

            for vtk_ele_idx, nid in enumerate( ele ):
                # node id is 1-based and the vtk index is 0-based
                vtk_ele.GetPointIds().SetId( vtk_ele_idx, nid - 1 )

            vtk_mesh.InsertNextCell( vtk_ele.GetCellType(), vtk_ele.GetPointIds() )

        return VTKMeshChecker( vtk_mesh )

    @classmethod
    def create_from_vtk_file( cls, file_name='tests.vtk', in_xml_format=True ):
        """
        Factory method: create checker from a given VTK mesh file

        :param file_name:
        :param in_xml_format:
        :return:
        """
        print( '--> Reading file \'{}\'...'.format( file_name ) )

        reader = vtk.vtkXMLUnstructuredGridReader() if in_xml_format else vtk.vtkUnstructuredGridReader( )
        reader.SetFileName( file_name )
        reader.Update()

        op = reader.GetOutput()

        vtk_mesh = vtk.vtkUnstructuredGrid()

        vtk_mesh.SetPoints( op.GetPoints() )

        # add the elements
        for ele_idx in range( op.GetNumberOfCells() ):
            vtk_ele = op.GetCell( ele_idx )
            vtk_mesh.InsertNextCell( vtk_ele.GetCellType(), vtk_ele.GetPointIds() )

        return VTKMeshChecker( vtk_mesh )

    def calculate_mesh_quality( self ):
        """
        The ranges were taken from
            http://www.vtk.org/Wiki/images/6/6b/VerdictManual-revA.pdf
        and
            http://csimsoft.com/help/trelishelp.htm
            -->  Mesh Generation / Mesh quality assessment / metrics for hexahedral elements

        Acceptable ranges are all *inclusive*
        - We only include metrics that are pertinent to a unit cube and not, for example, the volume of the
          whole mesh like 'RelativeSizeSquared' and others that use it.
        - We do not use the 'Oddy' metric because it gives values >0.5 for small cubes,
          which should be perfect, i.e. zero.

        Note 1: The definition of distortion varies slightly between the above two sources.
                The looser Verdict range (0.5, 1) is favoured over the Trelis range (0.6,1).

        Note 2: Use Trelis's 'Aspect Ratio' metric because the Verdict range of (1,1.3) applies to
                principal axes and it seems to be more appropriate to use edge length ratios (the Trelis defn).
                They both come from the same reference anyway (Taylor, 1989), which is unavailable.

        ***
        Code adapted from https://github.com/Kitware/VTK/blob/master/Filters/Verdict/Testing/Python/MeshQuality.py

        Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen.
        All rights reserved.
        See Copyright.txt or http://www.kitware.com/Copyright.htm for details.
        This software is distributed WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
        FITNESS FOR A PARTICULAR PURPOSE.
        See the above copyright notice for more information.
        ***
        """

        class HexMetricSummary(object):
            """ Metric summary for all VTK Hexahedra in mesh """

            def __init__( self, metric_name, valid_range ):
                self.metric_name = metric_name

                # store the given range...
                self.supplied_valid_range = valid_range

                #...but for checking we use a magic number
                lb = -1e99 if valid_range[0] is None else float( valid_range[0] )
                ub =  1e99 if valid_range[1] is None else float( valid_range[1] )
                self.valid_range = ( lb, ub )

                self.quality_values = dict( min=None, mean=None, max=None, sd=None )
                self.num_elements = None

                self.num_failures = 0

            def in_valid_range(self, value):
                """
                Check the range for a given value
                :param value: value of a metric
                :return: bool
                """
                if self.valid_range[0] <= value <= self.valid_range[1]:
                    return True

                # make sure the value is the correct sign if the lower bound is zero
                if value < self.valid_range[0 ] == 0.0:
                    return False

                # allow a margin around the extremes of the valid range...
                validity_margin = 1e-3

                if self.valid_range[1] + validity_margin < value or value < self.valid_range[0] - validity_margin:
                    return False

                # print( '--> {}: value is { }, valid range is {}'.format( self.metric_name, value, self.valid_range ) )

                return True

            def store_quality_values( self, mesh_quality_output ):
                """
                Get and store metrics from VTK

                :param mesh_quality_output: VTK mesh quality data
                :return:
                """

                cell_data = mesh_quality_output.GetCellData().GetArray( 'Quality' )
                field_data = mesh_quality_output.GetFieldData( ).GetArray( 'Mesh Hexahedron Quality' )

                # quality statistics are [ min, mean, max, variance, number of elements ]
                self.quality_values[ 'min' ]  = field_data.GetComponent( 0, 0 )
                self.quality_values[ 'mean' ] = field_data.GetComponent( 0, 1 )
                self.quality_values[ 'max' ]  = field_data.GetComponent( 0, 2 )
                self.quality_values[ 'sd' ]   = math.sqrt( math.fabs( field_data.GetComponent( 0, 3 ) ) )

                self.num_elements = int( field_data.GetComponent( 0, 4 ) )

                # for idx in range( qty.GetNumberOfTuples() ):
                for idx in range( self.num_elements ):
                    quality_value = cell_data.GetTuple1( idx )
                    if not self.in_valid_range( quality_value ):
                        self.num_failures += 1

            def __str__(self):
                format_string = 'Metric: {:20s} - {}. Valid range is {} - {}. ' + \
                                'Quality values: {} elements in {:5.3f} - {:5.3f}, mean = {:5.3f}, s.d. = {:5.3f}'

                return format_string.format(
                    self.metric_name,
                    'PASS' if self.num_failures == 0 else '# failure(s) = {} ({:.2f} %)'.format(
                        self.num_failures, 100 * float( self.num_failures ) / self.num_elements ),
                    self.supplied_valid_range[0],
                    self.supplied_valid_range[1],
                    self.num_elements,
                    self.quality_values[ 'min' ],
                    self.quality_values[ 'max' ],
                    self.quality_values[ 'mean' ],
                    self.quality_values[ 'sd' ]
                )

        metrics = [
            [ 'Condition',           (1,    8) ],
            [ 'Diagonal',            (0.65, 1) ],
            [ 'Distortion',          (0.5,  1) ],    # Note #1
            [ 'Jacobian',            (0,    None) ],
            [ 'MaxAspectFrobenius',  (1,    3) ],    # has it replaced 'condition'?
            [ 'MedAspectFrobenius',  (1,    3) ],
            [ 'MaxEdgeRatios',       (1,    4) ],    # Note #2
            [ 'ScaledJacobian',      (0.5,  1) ],
            [ 'Shape',               (0.3,  1) ],
            [ 'Shear',               (0.3,  1) ],
            [ 'Skew',                (0,    0.5) ],
            [ 'Stretch',             (0.25, 1) ],
            [ 'Taper',               (0,    0.5) ],
            [ 'Volume',              (0,    None ) ]
        ]

        print( '--> Checking mesh using {}'.format( vtk.vtkVersion.GetVTKSourceVersion() ) )
        print( '--> *ONLY* hexahedral element quality is calculated by vtk - other element types are ignored!' )

        # Set up the mesh quality object
        iq = vtk.vtkMeshQuality()
        iq.SetInputData( self.vtk_mesh )

        quality_results = list()

        for metric in metrics:
            hex_metric = HexMetricSummary( metric[0], metric[1] )

            cmd = 'iq.SetHexQualityMeasureTo{}()'.format( hex_metric.metric_name )
            eval( cmd )

            iq.Update()

            hex_metric.store_quality_values( iq.GetOutput() )

            quality_results.append( str( hex_metric ) )

        return quality_results

    def output_vtk_file(self, file_name='tests.vtk', in_xml_format=True ):
        """
        Write the mesh to a VTK file

        :param file_name:
        :param in_xml_format:
        :return:
        """

        gw = vtk.vtkXMLUnstructuredGridWriter() if in_xml_format else vtk.vtkUnstructuredGridWriter()

        if in_xml_format:
            # XML data written in binary format by default
            gw.SetDataModeToAscii()

        gw.SetFileName( file_name )
        gw.SetInputData( self.vtk_mesh )

        rtc = gw.Write()

        if rtc == 1:
            print( '--> Wrote mesh to {}'.format( file_name ) )
        else:
            print( '--> Failed to write mesh to file')

        # Exodus format
        # exow = vtk.vtkExodusIIWriter()
        # exow.SetFileName( 'tests.exo' )
        # exow.SetInputData( mesh )
        #
        # rtc = exow.Write()

        return

    def view_mesh(self):
        """
        Use VTK to visualise the mesh

        Note: colours are given in R, G, B format with each number in [0, 1]

        :return:
        """

        mapper = vtk.vtkDataSetMapper()
        mapper.SetInputData( self.vtk_mesh )

        actor = vtk.vtkActor()
        actor.SetMapper( mapper )

        # set element colour
        actor.GetProperty().SetColor(0.7, 0.9, 0.8)

        # set edge colour
        actor.GetProperty().SetEdgeColor( 0.1, 0.1, 0.1 )
        actor.GetProperty().EdgeVisibilityOn()

        renderer = vtk.vtkRenderer()
        renderer.AddActor( actor )
        renderer.SetBackground(0.2, 0.3, 0.4)

        render_window = vtk.vtkRenderWindow()
        render_window.AddRenderer( renderer )
        render_window.Render()

        render_window_interactor = vtk.vtkRenderWindowInteractor()
        render_window_interactor.SetRenderWindow( render_window )
        render_window_interactor.Start()

        return


def make_simple_vtk_mesh( ):
    """ Make a simple VTK mesh """

    point_array = list()

    def from_xyz( x, y, z ):
        """ helper function to construct a point """
        return gp.Point( x, y, z, allocate_id=True )

    point_array.append( from_xyz( 0.0, 0.2, 0.0 ) )
    point_array.append( from_xyz( 1.1, 0.1, 0.0 ) )
    point_array.append( from_xyz( 1.0, 1.2, 0.0 ) )
    point_array.append( from_xyz( 0.0, 1.3, 0.0 ) )

    # in z=1
    point_array.append( from_xyz( 0.0, 0.0, 1.0 ) )
    point_array.append( from_xyz( 1.2, 0.0, 1.0 ) )
    point_array.append( from_xyz( 1.0, 1.1, 1.0 ) )
    point_array.append( from_xyz( 0.0, 1.0, 1.0 ) )

    # in z=-1
    point_array.append( from_xyz( 0.0, 0.0, -1.0 ) )
    point_array.append( from_xyz( 0.9, 0.0, -1.0 ) )
    point_array.append( from_xyz( 1.0, 1.0, -1.0 ) )
    point_array.append( from_xyz( 0.0, 1.0, -1.0 ) )

    # for the pentahedron
    point_array.append( from_xyz( 0.0, 1.0, 2.0 ) )
    point_array.append( from_xyz( 0.0, 0.0, 2.0 ) )


    elements = ( range(1, 9), (9, 10, 11, 12, 1, 2, 3, 4), (8, 7, 13, 5, 6, 14) )

    return VTKMeshChecker.create_from_points_and_elements( point_array, elements )

def make_quarter_disc_vtk_hex_mesh():
    """
    Makes a quarter-circle mesh
    """

    # Create the points

    from numpy import sin, pi, sqrt, linspace

    num_r = 10

    dr = 0.2
    dz = 0.2

    sqrt2 = sqrt( 2.0 )

    point_array = list()

    # central point
    point_array.append( gp.Point( 0, 0, 0 ) )

    for ri in range( num_r ):
        num_pts = 2 * ri + 3

        sint = [ sin(t) for t in linspace(0.0, pi/2, num_pts) ]
        cost = reversed( sint )

        for st, ct in zip( sint, cost ):
            x = (ri + 1) * dr * ( ct + min( 1.0, sqrt2 * ct ) ) / 2
            y = (ri + 1) * dr * ( st + min( 1.0, sqrt2 * st ) ) / 2

            point_array.append( gp.Point( x, y, 0 ) )

    num_pts = len( point_array )

    dp = ( 0, 0, dz )
    point_array = point_array + [ gp.translate_pt_by_xyz( pt, dp ) for pt in point_array ]

    for pt in point_array:
        pt.allocate_id()

    # Create elements

    elements = list()

    for ri in range( num_r ):
        riri = ri*ri
        tworip1 = 2*ri + 1
        range_1 = range( riri, riri + tworip1 )
        range_2 = range( range_1[-1] + 1, range_1[-1] + tworip1 + 3 )

        for idx, ele_num in enumerate( range_1 ):
            if 2 * idx + 1 == len( range_1 ):
                array_indices = [ ele_num, ] + range_2[ idx:idx+3 ]
            elif 2 * idx < len( range_1 ):
                array_indices = [ ele_num, range_2[idx], range_2[idx+1], range_1[idx+1] ]
            else:
                array_indices = [ ele_num - 1, range_2[idx+1], range_2[idx+2], range_1[idx] ]

            array_indices += [ ai + num_pts for ai in array_indices ]

            elements.append( [ point_array[arr_idx].id for arr_idx in array_indices] )

    return VTKMeshChecker.create_from_points_and_elements( point_array, elements )

def main():
    """ Main method """

    option = 1
    mesh_checker = None

    if option == 1:
        mesh_checker = make_simple_vtk_mesh( )

    if option == 2:
        mesh_checker = make_quarter_disc_vtk_hex_mesh()

    if option == 3:
        mesh_checker = VTKMeshChecker.create_from_vtk_file( 'tests.vtk' )

    # run the mesh checker
    quality_results = mesh_checker.calculate_mesh_quality( )
    for qr in quality_results:
        print( '--> {}'.format( qr ) )

    mesh_checker.output_vtk_file()
    mesh_checker.view_mesh()

if __name__ == '__main__':
    main()
