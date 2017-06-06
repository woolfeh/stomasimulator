"""
Write the 'feb' file for the FEBio simulation


Example feb file - not an exhaustive list of features

<?xml version="1.0" encoding="ISO-8859-1"?>
<febio_spec version="2.0">
  <Module type="solid"/>
  <Control>
    <time_steps>50</time_steps>
    <step_size>0.1</step_size>
    <max_refs>15</max_refs>
    <max_ups>10</max_ups>
    <dtol>0.001</dtol>
    <etol>0.01</etol>
    <rtol>0</rtol>
    <lstol>0.9</lstol>
    <time_stepper>
      <dtmin>0.01</dtmin>
      <dtmax>0.05</dtmax>
      <max_retries>10</max_retries>
      <opt_iter>10</opt_iter>
      <aggressiveness>1</aggressiveness>
    </time_stepper>
    <analysis type="static"/>
  </Control>
  <Globals>
    <Constants>
      <T>0</T>
      <R>0</R>
      <Fc>0</Fc>
    </Constants>
  </Globals>
  <Material>
    <material id="1" name="MaterialPosX" type="trans iso Mooney-Rivlin">
      <density>1</density>
      <c1>5</c1>
      <c2>5</c2>
      <k>10000</k>
      <c3>0</c3>
      <c4>0</c4>
      <c5>0</c5>
      <lam_max>0</lam_max>
      <fiber type="local">     3,     4</fiber>
    </material>
    <material id="2" name="MaterialPosXPole" type="trans iso Mooney-Rivlin">
      <density>1</density>
      <c1>5</c1>
      <c2>5</c2>
      <k>10000</k>
      <c3>0</c3>
      <c4>0</c4>
      <c5>0</c5>
      <lam_max>0</lam_max>
      <fiber type="cylindrical">
        <center>18.5625,0,0</center>
        <axis>0,1,0</axis>
        <vector>0,0,1</vector>
      </fiber>
    </material>
  </Material>
  <Geometry>
    <Nodes>
      <node id="1">  1.1210000e+01,  0.0000000e+00,  0.0000000e+00</node>
      ...
      <node id="40000"> -1.3617659e+01,  1.5269246e-01, -3.5319526e-01</node>
    </Nodes>
    <Elements type="hex8" mat="1" elset="Part1">
      <elem id="1"> 14823,  2402,    87,  2414, 34419, 15747,  2514, 15411</elem>
      ...
      <elem id="40000"> 43511,  4298, 25058, 22655,   143,  3982</elem>
    </Elements>
  </Geometry>
  <Boundary>
    <fix bc="x">
      <node id="75"/>
      ...
      <node id="43553"/>
    </fix>
  </Boundary>
  <Loads>
    <surface_load type="pressure">
      <pressure lc="1">1</pressure>
      <linear>0</linear>
      <surface>
        <quad4 id="1">   565,   566,  5832,  5811</quad4>
        ...
        <quad4 id="1056">  6818,   564,    15,   576</quad4>
        <tri3 id="1057">    14,   565,  5811</tri3>
        ...
        <tri3 id="1144">    14,  6587,   565</tri3>
      </surface>
    </surface_load>
  </Loads>
  <Contact>
    <contact type="rigid_wall">
      <laugon>1</laugon>
      <tolerance>0.1</tolerance>
      <penalty>10</penalty>
      <plane lc="2">0,1,0,0</plane>
      <surface>
        <quad4 id="1">   385,    10,   613,  6819</quad4>
        ...
        <quad4 id="7392">  4290, 24950,  3885,   139</quad4>
      </surface>
    </contact>
  </Contact>
  <LoadData>
    <loadcurve id="1" type="linear">
      <point>0,0</point>
      <point>5,5</point>
    </loadcurve>
    <loadcurve id="2" type="linear">
      <point>0,0</point>
      <point>5,0</point>
    </loadcurve>
  </LoadData>
  <Output>
    <plotfile type="febio">
      <var type="displacement"/>
      <var type="stress"/>
    </plotfile>
  </Output>
</febio_spec>

"""

from __future__ import print_function

from lxml import etree

import stomasimulator.mesh.stoma_mesh as sm
import stomasimulator.stomata.stoma_config as sc

from . import feb_file_classes as ffc
from . import feb_material as febmat

# TODO: Use a separate file for the material definition

# This can be accomplished by changing the 'Material' section in the main file to:
#
# <Material from="name-of-material-definition-file.feb"/>
#
# and where an example 'name-of-material-definition-file.feb' could be:
#
# <?xml version='1.0' encoding='utf-8'?>
# <febio_spec version="2.0">
#   <Material>
#     <material id="1" name="Material1" type="trans iso Veronda-Westmann">
#       <density>1</density>
#       <c1>2.0</c1>
#       <c2>2.0</c2>
#       <k>10000.0</k>
#       <c3>0.0</c3>
#       <c4>0.0</c4>
#       <c5>1000.0</c5>
#       <lam_max>1.0</lam_max>
#     </material>
#   </Material>
# </febio_spec>


def write_control_section( root, timestep_params ):
    """
      <Control>
        <time_steps>50</time_steps>
        <step_size>0.1</step_size>
        <max_refs>15</max_refs>
        <max_ups>10</max_ups>
        <dtol>0.001</dtol>
        <etol>0.01</etol>
        <rtol>0</rtol>
        <lstol>0.9</lstol>
        <time_stepper>
          <dtmin>0.01</dtmin>
          <dtmax>0.05</dtmax>
          <max_retries>10</max_retries>
          <opt_iter>10</opt_iter>
          <aggressiveness>1</aggressiveness>
        </time_stepper>
        <analysis type="static"/>
      </Control>
    """

    control_element = etree.SubElement( root, 'Control' )

    etree.SubElement( control_element, 'time_steps' ).text = '{}'.format( timestep_params.num_steps )
    etree.SubElement( control_element, 'step_size' ).text = '{}'.format( timestep_params.dt )

    # use BFGS
    # etree.SubElement( control_element, 'max_refs'   ).text = '15'
    # etree.SubElement( control_element, 'max_ups'    ).text = '10'

    # use full-Newton method (more stiffness matrix reformations)
    etree.SubElement( control_element, 'max_refs' ).text = '50'
    etree.SubElement( control_element, 'max_ups' ).text = '0'

    etree.SubElement( control_element, 'dtol' ).text = '0.001'
    etree.SubElement( control_element, 'etol' ).text = '0.01'
    etree.SubElement( control_element, 'rtol' ).text = '0'
    etree.SubElement( control_element, 'lstol' ).text = '0.9'

    ts_element = etree.SubElement( control_element, 'time_stepper' )
    etree.SubElement( ts_element, 'dtmin' ).text = '{}'.format( timestep_params.min_dt )
    etree.SubElement( ts_element, 'dtmax' ).text = '{}'.format( timestep_params.max_dt )
    etree.SubElement( ts_element, 'max_retries' ).text = '15'
    etree.SubElement( ts_element, 'opt_iter' ).text = '10'
    etree.SubElement( ts_element, 'aggressiveness' ).text = '1'

    etree.SubElement( control_element, 'analysis', type='static' )


def write_globals_section( root ):
    """
    <Globals>
        <Constants>
            <T>0</T>
            <R>0</R>
            <Fc>0</Fc>
        </Constants>
    </Globals>

    :return:
    """

    globals_element = etree.SubElement( root, 'Globals' )
    constants_element = etree.SubElement( globals_element, 'Constants' )

    for const in ('T', 'R', 'Fc'):
        etree.SubElement( constants_element, const ).text = '0'


def write_material_section( root, material_model ):
    """
    <Material>
      <material id="1" name="Material1" type="trans iso Mooney-Rivlin">
        <density>1</density>
        <c1>5</c1>
        <c2>5</c2>
        <k>10000</k>
        <c3>0</c3>
        <c4>0</c4>
        <c5>0</c5>
        <lam_max>0</lam_max>
      </material>
    </Material>

    The fibre vector can be added to the material definition using,
      <fiber type="local">3,4</fiber>
    for example, but we prefer the per-element definition set in ElementData

    :param root: root XML node
    :param material_model:
    :type material_model: febmat.MaterialModel
    :return:
    """

    def add_material( mat_id ):
        """ Add the material to the XML """

        material_element = etree.SubElement( materials_element, 'material',
                                             id='{}'.format( mat_id ),
                                             name='Material{}'.format( mat_id ),
                                             type=material_model.febio_name )

        etree.SubElement( material_element, 'density' ).text = '1'

        matrix = material_model.isotropic_matrix

        etree.SubElement( material_element, 'c1' ).text = '{}'.format( matrix.c1 )
        etree.SubElement( material_element, 'c2' ).text = '{}'.format( matrix.c2 )
        etree.SubElement( material_element, 'k' ).text = '{}'.format( matrix.bulk_modulus( ) )

        if not material_model.is_isotropic:
            fibres = material_model.fibres

            etree.SubElement( material_element, 'c3' ).text = '{}'.format( fibres.c3 )
            etree.SubElement( material_element, 'c4' ).text = '{}'.format( fibres.c4 )
            etree.SubElement( material_element, 'c5' ).text = '{}'.format( fibres.c5 )
            etree.SubElement( material_element, 'lam_max' ).text = '{}'.format( fibres.lm )

        return

    # Create the main materials elements
    materials_element = etree.SubElement( root, 'Material' )

    # Add materials to the node
    add_material( 1 )


def write_geometry_section( root, mesh ):
    """ Write the geometry to the XML

    <Geometry>
        <Nodes>...</Nodes>
        <Elements>...</Elements>
        <Elements>...</Elements>
        ...

        <ElementData>...</ElementData>
    </Geometry>
    """

    def write_nodes_section( ):
        """
        <Nodes>
          <node id="1">  1.1210000e+01,  0.0000000e+00,  0.0000000e+00</node>
          ...
          <node id="43742"> -1.3617659e+01,  1.5269246e-01, -3.5319526e-01</node>
        </Nodes>
        """

        nodes_element = etree.SubElement( geometry_element, 'Nodes' )

        # the nodes must be added with a sequential id - this assumes 'nodes_map' fulfils that criterion
        for node_id, node in mesh.nodes_map.items( ):
            node_str = '{}, {}, {}'.format( node.x, node.y, node.z )
            etree.SubElement( nodes_element, 'node', id='{}'.format( node_id ) ).text = node_str

    def write_elements_section( ):
        """
        Add each of the elements - partitioned by element type, material and part

        <Elements type="hex8" mat="1" elset="Part1">
            <elem id="1"> 14823,  2402,    87,  2414, 34419, 15747,  2514, 15411</elem>
            <elem id="2"> 34419, 15747,  2514, 15411, 34420, 15748,  2515, 15432</elem>
            ...
        </Elements>
        <Elements type="penta6" mat="2" elset="Part2">
            ...
        </Elements>
        """

        idx = 1
        for ele_type, ele_list in mesh.get_elements_by_type( ).items( ):
            elements_node = etree.SubElement( geometry_element, 'Elements', type=ele_type, mat='1', elset='Part1' )

            for ele in ele_list:
                elem_node = etree.SubElement( elements_node, 'elem', id='{}'.format( idx ) )
                elem_node.text = ele.formatted_node_list

                idx += 1

    def write_element_data_section( ):
        """
        Add the fibre direction to
        o  This XML section must be after the 'Elements' section (otherwise PreView crashes)
        o  Use the index of addition for the id rather than the element's id because it just needs to be sequential...
           ...so iterate the same way as in write_elements_section

        <ElementData>
          <element id="1"><fiber>0.1,0.2,0.3</fiber>
          ...
        </ElementData>
        """

        element_data_se = etree.SubElement( geometry_element, 'ElementData' )

        idx = 1
        for _, ele_list in mesh.get_elements_by_type( ).items( ):
            for ele in ele_list:
                fv = ele.fibre_vector
                if fv is not None:
                    ele_se = etree.SubElement( element_data_se, 'element', id='{}'.format( idx ) )
                    etree.SubElement( ele_se, 'fiber' ).text = '{},{},{}'.format( fv.x, fv.y, fv.z )

                idx += 1

        return

    geometry_element = etree.SubElement( root, 'Geometry' )

    write_nodes_section( )
    write_elements_section( )
    write_element_data_section()

    return


def write_bc_section( root, mesh ):
    """
    Write the boundary conditions

    <Boundary>
      <fix bc="x"><node id="75"/>...</fix>
      <fix bc="y"><node id="43"/>...</fix>
      <fix bc="z"><node id="12"/>...</fix>
    </Boundary>

    :param root:
    :param mesh:
    :type mesh: sm.SimpleMesh
    :return:
    """

    boundary_element = etree.SubElement( root, 'Boundary' )

    # Add the nodes for each boundary condition
    #
    for fixed_dof, node_list in mesh.bcs.items( ):
        if len( node_list ) == 0:
            continue

        bc_element = etree.SubElement( boundary_element, 'fix', bc=fixed_dof )

        for node_id in node_list:
            etree.SubElement( bc_element, 'node', id='{}'.format( node_id ) )


def write_loads_section( root, mesh, load_curves ):
    """
    Write the pressure loads

    <Loads>
      <surface_load type="pressure">
        <pressure lc="1">1</pressure>
        <linear>0</linear>
        <surface>
          <quad4 id="1">   565,   566,  5832,  5811</quad4>
          ...
          <tri3 id="9327"> 22654,   143, 22655</tri3>
          ...
        </surface>
      </surface_load>
    </Loads>

    :param root:
    :param mesh:
    :param load_curves:
    :param timestep_params:
    :return:
    """

    def add_pressure_facets( end_pressure, pressure_facets ):
        """
        Add pressure facets to the XML document
        :param end_pressure:
        :param pressure_facets:
        :return:
        """

        surface_load_element = etree.SubElement( loads_node, 'surface_load', type='pressure' )

        # create the load curve (1-to-1 scaling between simulation time and stoma pressure)
        lc = ffc.LoadCurve( )
        lc.add_point( 0, 0 )
        lc.add_point( end_pressure, end_pressure )
        load_curves.append( lc )

        # scaling factor
        etree.SubElement( surface_load_element, 'pressure', lc='{}'.format( lc.load_curve_id ) ).text = '1'
        etree.SubElement( surface_load_element, 'linear' ).text = '0'

        surface_element = etree.SubElement( surface_load_element, 'surface' )

        for idx, facet in enumerate( pressure_facets ):
            # id for the facet is 1-based
            add_facet( parent_element=surface_element, facet=facet, idx=idx + 1 )


    internal_facets = sum( len( f ) for f in mesh.pressure_facets )

    stoma_pressure = get_max_stoma_pressure( mesh )
    dorsal_pressure = get_max_dorsal_pressure( mesh )
    num_gcs = mesh.mesh_config.stoma_cfg.num_gcs

    if internal_facets > 0 or dorsal_pressure > 0.0:
        loads_node = etree.SubElement( root, 'Loads' )

    if internal_facets > 0:
        for gc_id in range( num_gcs ):
            add_pressure_facets( end_pressure=stoma_pressure,
                                 pressure_facets=mesh.pressure_facets[ gc_id ] )

    if get_max_dorsal_pressure( mesh ) > 0.0:
        for gc_id in range( num_gcs ):
            add_pressure_facets( end_pressure=dorsal_pressure,
                                 pressure_facets=mesh.dorsal_wall_facets[ gc_id ] )

    return


def write_contact_section( root, mesh, load_curves ):
    """
    <Contact>
      <contact type="rigid_wall">
        <laugon>1</laugon>
        <tolerance>0.1</tolerance>
        <penalty>10</penalty>
        <plane lc="2">0,1,0,0</plane>
        <surface>
          <quad4 id="1">   385,    10,   613,  6819</quad4>
          ...
          <quad4 id="7392">  4290, 24950,  3885,   139</quad4>
        </surface>
      </contact>
    </Contact>

    :param root:
    :param mesh:
    :param load_curves:
    :return:
    """

    def add_facets_to_surface( surface_xml_ele, facet_list ):
        """ Add the facets to a surface """
        local_idx = 1
        for facet in facet_list:
            add_facet( surface_xml_ele, facet, local_idx )
            local_idx += 1

    def add_sliding_contact_def():
        """ add sliding interface """

        def_ele = etree.SubElement( contact_element, 'contact', type='facet-to-facet sliding' )

        etree.SubElement( def_ele, 'laugon' ).text = '0'
        etree.SubElement( def_ele, 'tolerance' ).text = '0.01'
        etree.SubElement( def_ele, 'penalty' ).text = '10'
        etree.SubElement( def_ele, 'two_pass' ).text = '0'
        etree.SubElement( def_ele, 'auto_penalty' ).text = '1'
        etree.SubElement( def_ele, 'search_tol' ).text = '0.01'
        etree.SubElement( def_ele, 'minaug' ).text = '0'
        etree.SubElement( def_ele, 'maxaug' ).text = '100'
        etree.SubElement( def_ele, 'gaptol' ).text = '0'

        return def_ele

    total_facets = sum( len( f ) for f in mesh.ventral_wall_facets.values( ) )

    num_gcs = mesh.mesh_config.stoma_cfg.num_gcs

    if total_facets > 0:
        # add the contact condition
        # ** only really required when the cell wall matrix is isotropic **
        contact_element = etree.SubElement( root, 'Contact' )

        if num_gcs == 1:
            # 1 GC so need contact with rigid wall to stop negative apertures

            # load curve that stops the rigid wall from moving
            lc = ffc.LoadCurve( )
            lc.add_point( 0, 0 )
            lc.add_point( get_max_stoma_pressure( mesh ), 0 )
            load_curves.append( lc )

            # add rigid wall (at y=0, x and z are free)
            #
            wall_contact_element = etree.SubElement( contact_element, 'contact', type='rigid_wall' )
            etree.SubElement( wall_contact_element, 'laugon' ).text = '1'
            etree.SubElement( wall_contact_element, 'tolerance' ).text = '0.1'
            etree.SubElement( wall_contact_element, 'penalty' ).text = '10'
            etree.SubElement( wall_contact_element, 'plane', lc='{}'.format( lc.load_curve_id ) ).text = '0,1,0,0'

            wall_contact_surface_element = etree.SubElement( wall_contact_element, 'surface' )

            add_facets_to_surface( surface_xml_ele=wall_contact_surface_element,
                                   facet_list=mesh.ventral_wall_facets[ 0 ] )

        if num_gcs == 2:
            # 2 GCs - need to make sure they stay apart

            # add sliding interface
            #
            f2f_sliding_element = add_sliding_contact_def()

            # master and slave will be the same, i.e. the ventral walls of both GCs
            for gc_id, surface_type in enumerate( ('master', 'slave') ):
                wall_contact_surface_element = etree.SubElement( f2f_sliding_element, 'surface', type=surface_type )

                add_facets_to_surface( surface_xml_ele=wall_contact_surface_element,
                                       facet_list=mesh.ventral_wall_facets[ gc_id ] )

    if len( mesh.epidermis_facets ) > 0:
        # the epidermal wall is present so add contact for that too

        f2f_sliding_element = add_sliding_contact_def()

        # master will be the epidermal wall
        contact_surface_element = etree.SubElement( f2f_sliding_element, 'surface', type='master' )
        add_facets_to_surface( contact_surface_element, mesh.epidermis_facets )

        # slave will be the dorsal walls
        contact_surface_element = etree.SubElement( f2f_sliding_element, 'surface', type='slave' )

        if num_gcs == 1:
            add_facets_to_surface( surface_xml_ele=contact_surface_element,
                                   facet_list=mesh.dorsal_wall_facets[ 0 ] )
        else:
            add_facets_to_surface( surface_xml_ele=contact_surface_element,
                                   facet_list=mesh.dorsal_wall_facets[ 0 ] + mesh.dorsal_wall_facets[ 1 ] )

    return


def write_load_data_section( root, load_curves ):
    """
    <LoadData>
      <loadcurve id="1" type="linear">
        <point>0,0</point>
        <point>5,5</point>
      </loadcurve>
      <loadcurve id="2" type="linear">
        <point>0,0</point>
        <point>5,0</point>
      </loadcurve>
    </LoadData>

    :param root:
    :param load_curves:
    :return:
    """

    load_data_element = etree.SubElement( root, 'LoadData' )

    for load_curve in load_curves:
        # type in {'step', 'linear', 'smooth'}
        ele = etree.SubElement( load_data_element,
                                'loadcurve',
                                id='{}'.format( load_curve.load_curve_id ),
                                type='linear' )

        for pair in load_curve.values:
            etree.SubElement( ele, 'point' ).text = '{}, {}'.format( pair[ 0 ], pair[ 1 ] )


def write_output_section( root ):
    """
    <Output>
      <plotfile type="febio">
        <var type="displacement"/>
        <var type="stress"/>
      </plotfile>
    </Output>

    :param root:
    :return:
    """

    output_element = etree.SubElement( root, 'Output' )
    plotfile_element = etree.SubElement( output_element, 'plotfile', type='febio' )

    for output_type in ( 'displacement', 'stress', 'fiber vector', 'relative volume' ):
        etree.SubElement( plotfile_element, 'var', type=output_type )


def add_facet( parent_element, facet, idx ):
    """ Add a surface facet to the XML

    :param parent_element:
    :param facet:
    :param idx:
    :return:
    """

    if facet is None:
        raise ValueError( 'The facet is bobbins - it is None' )

    if len( facet ) not in (3, 4):
        raise ValueError( 'The facet is bobbins - its length should be 3 or 4 but it is {}'.format( len( facet ) ) )

    # Formatting is FEBio-specific so handle it here (could do it in the mesh)

    if len( facet ) == 3:
        facet_type_name = 'tri3'
        format_string = '{}, {}, {}'
    else:
        facet_type_name = 'quad4'
        format_string = '{}, {}, {}, {}'

    facet_element = etree.SubElement( parent_element, facet_type_name, id='{}'.format( idx ) )
    facet_element.text = format_string.format( *facet )

    return


def get_max_stoma_pressure( mesh ):
    """
    :param mesh:
    :return:
    """
    return mesh.mesh_config.stoma_cfg.pressure.stoma_max

def get_max_dorsal_pressure( mesh ):
    """
    :param mesh:
    :return:
    """
    return mesh.mesh_config.stoma_cfg.pressure.epidermis_max

def create_feb_xml( mesh ):
    """ Create the XML for the 'feb' file

    :param mesh:
    :return
    """

    ffc.LoadCurve.reset_load_curve_id( )

    load_curves = [ ]

    timestep_params = ffc.get_timestep_control( end_time=get_max_stoma_pressure( mesh ) )

    # create the root element (must have the spec version - we will use v2)
    root = etree.Element( 'febio_spec', version='2.0' )

    # add the type of analysis
    etree.SubElement( root, 'Module', type='solid' )

    write_control_section( root, timestep_params )
    write_globals_section( root )
    write_material_section( root, mesh.mesh_config.stoma_cfg.material_model )
    write_geometry_section( root, mesh )
    write_bc_section( root, mesh )
    write_loads_section( root, mesh, load_curves )
    write_contact_section( root, mesh, load_curves )
    write_load_data_section( root, load_curves )
    write_output_section( root )

    xml_str = etree.tostring( element_or_tree=root,
                              encoding='utf-8',
                              xml_declaration=True,
                              pretty_print=True )

    return xml_str


def output_febio_xml( stoma_cfg, filename, perform_mesh_checks=True ):
    """ Create and mesh and output it to a 'feb' file

    :param stoma_cfg:
    :type stoma_cfg: sc.StomaConfig
    :param filename:
    :param perform_mesh_checks:
    :return:
    """

    stoma_cfg.print_detail( )

    mesh = sm.build_stoma_mesh( stoma_cfg, perform_mesh_checks=perform_mesh_checks )

    xml_str = create_feb_xml( mesh )

    with open( filename, 'w' ) as fh:
        fh.write( xml_str )

    print( '--> FEBio file written to {}'.format( filename ) )

    return filename


if __name__ == '__main__':
    raise NotImplementedError( 'Access to this functionality is via the top-level module' )
