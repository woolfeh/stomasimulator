""" Empty module that acts as a stub for an include in mesh_check """

from __future__ import print_function

# pylint: disable=missing-docstring,invalid-name

def log_message( msg ):
    print( 'vtk_empty: {}'.format( msg ) )


class Dummy(object):
    """ A catch-all object for the null interface """

    @staticmethod
    def GetVTKSourceVersion( ):
        return ''

    def SetId(self, i, j):
        pass

    def SetColor( self, aa, bb, cc ):
        pass

    def SetEdgeColor( self, aa, bb, cc, ):
        pass

    def EdgeVisibilityOn( self ):
        pass

DUMMY = Dummy()

vtkVersion = DUMMY

class vtkPoints( object ):
    @staticmethod
    def InsertNextPoint( _ ):
        log_message( 'Ignoring InsertNextPoint' )


class vtkUnstructuredGrid( object ):
    @staticmethod
    def SetPoints( _ ):
        log_message( 'Ignoring call to SetPoints' )

    @staticmethod
    def InsertNextCell( aa, bb ):
        log_message( 'Ignoring call to InsertNextCell' )


class vtkUnstructuredGridReader( object ):
    pass


class vtkUnstructuredGridWriter( object ):
    pass


class vtkHexahedron(object):
    @staticmethod
    def GetPointIds( ):
        return DUMMY

    @staticmethod
    def GetCellType( ):
        return ''


class vtkWedge(vtkHexahedron):
    pass


class vtkDataSetMapper(object):
    def SetInputData( self, a ):
        pass


class vtkActor(object):
    def SetMapper( self, a ):
        pass

    @staticmethod
    def GetProperty( ):
        return DUMMY

class vtkRenderer(object):
    def AddActor( self, a ):
        pass

    def SetBackground( self, a, b, c ):
        pass


class vtkRenderWindow(object):
    def AddRenderer( self, a ):
        pass

    def Render( self ):
        pass


class vtkRenderWindowInteractor(object):
    def SetRenderWindow( self, a ):
        pass

    def Start(self):
        pass


class vtkXMLUnstructuredGridReader(object):
    def SetFileName(self, a):
        pass

    def Update(self):
        pass

    def GetOutput(self):
        pass


class vtkXMLUnstructuredGridWriter(object):
    def SetDataModeToAscii(self):
        pass

    def SetFileName( self, a ):
        pass

    def SetInputData( self, a ):
        pass

    @staticmethod
    def Write( ):
        return 0


class vtkMeshQuality(object):
    def SetInputData(self,a):
        pass

    def Update(self):
        pass

    @staticmethod
    def GetOutput( ):
        return ''


if __name__ == '__main__':
    pass
