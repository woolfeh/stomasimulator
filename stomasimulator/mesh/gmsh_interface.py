""" Write mesh out in Gmsh format """

from __future__ import print_function

def write_to_file(mesh, file_name):
    """
    Write the given mesh to a Gmsh file
    :param mesh: A SimpleMesh object
    :param file_name: name of the Gmsh file
    :return:
    """

    with open( file_name, 'w' ) as fh:
        for pt in mesh.nodes_map.values():
            if pt is None:
                continue

            # Format for a point is: Point(pt_id) = { x, y, z };
            fh.write( 'Point({}) = {{ {}, {}, {} }};\n'.format( pt.id, pt.x, pt.y, pt.z ) )

        fh.write( '\n' )

        # Format for an edge is: Line(edge_id) = { pt_id_1, pt_id_2 }
        for edge_idx, edge in enumerate( mesh.edges ):
            fh.write( 'Line({}) = {{ {}, {} }};\n'.format( edge_idx + 1, edge[0], edge[1] ) )

    print( '--> Written Gmsh file {}'.format( file_name ) )

if __name__ == '__main__':
    pass
