
import matplotlib as mpl
import matplotlib.pyplot as pplt

import stomasimulator.geom.ellipse as el
import stomasimulator.geom.line as line
import stomasimulator.geom.point as pt
import stomasimulator.geom.geom_utils as g

from . import stoma_mesh as msh


def plot_equatorial_points( cfg ):

    nnn = cfg.mesh_config.num_slices

    ellipse_d = el.Ellipse( cfg.a_d, cfg.b_d )
    ellipse_v = el.Ellipse( cfg.a_v, cfg.b_v )

    # Calculate the points on the ventral and dorsal walls
    xy_v = el.calculate_pts_for_equi_arc_spaced_t( ellipse_v, nnn )
    xy_d = el.calculate_pts_for_equi_arc_spaced_t( ellipse_d, nnn )

    # get the 'real' PM points s.t. the wall thickness is conserved
    xy_v_pm_2 = el.calculate_polyline_offset_from_ellipse( ellipse_v,  cfg.wall_thickness.ventral )
    xy_d_pm_2 = el.calculate_polyline_offset_from_ellipse( ellipse_d, -cfg.wall_thickness.dorsal  )

    # calculate equal-polar angle spaced points
    xy_th1 = el.calculate_pts_for_equi_spaced_t( ellipse_d, nnn )
    xy_th2 = el.calculate_pts_for_equi_spaced_t( ellipse_v, nnn )

    plt, ax = el.plot_initialise( )

    # plot the lines between the equi-arc spaced points
    for pts in zip( xy_d, xy_v ):
        plt.plot( (pts[0].x, pts[1].x), (pts[0].y, pts[1].y), 'r' )

    # plot the lines between the equi-polar angle spaced points
    for pts in zip( xy_th1, xy_th2 ):
        plt.plot( (pts[0].x, pts[1].x), (pts[0].y, pts[1].y), 'g' )

    for pts, pt_spec in zip( ( xy_v, xy_d, xy_v_pm_2, xy_d_pm_2 ),
                             ( '.r', '.r', '.r',      '.r' ) ):
        el.plot_add_pts( plt, pts, pt_spec )

    for a, b in zip( ( cfg.a_mid, cfg.a_d, cfg.a_v ),
                     ( cfg.b_mid, cfg.b_d, cfg.b_v ) ):
        ellipse = el.Ellipse( a, b )
        some_pts = el.calculate_pts_for_equi_spaced_t( ellipse, 100 )
        el.plot_add_pts( plt, some_pts, '-k' )

    for pt_pair in zip( xy_v, xy_d ):
        l = line.Line( pt_pair[0], pt_pair[1] )

        # find ventral PM points from v-d line intersection with offset points
        poi_v = g.find_point_of_intersection( l, xy_v_pm_2 )

        if poi_v is not None:
            plt.plot( poi_v.x, poi_v.y, 'ob' )

        # find dorsal PM points from v-d line intersection with offset points
        poi_d = g.find_point_of_intersection( l, xy_d_pm_2 )

        if poi_d is not None:
            plt.plot( poi_d.x, poi_d.y, 'ob' )

        if poi_d is not None and poi_v is not None:
            mid_pt = 0.5 * ( poi_d + poi_v )
            plt.plot( mid_pt.x, mid_pt.y, 'ob' )

    # draw a line based on the most troublesome point
    mtp = pt.Point( cfg.a_v + cfg.wall_thickness.ventral, cfg.wall_thickness.polar, 0.0 )
    mtp_n = mtp.unit()

    plt.plot( mtp.x, mtp.y, 'ok')
    plt.plot( (0.0, cfg.a_d * mtp_n.x), (0.0, cfg.a_d * mtp_n.y), 'k')

    plt.show()

    return


def plot_mesh( mesh, add_labels=False ):
    # need to load this - required by the subplot projection
    import mpl_toolkits.mplot3d.axes3d as ax3

    mpl.rcParams['legend.fontsize'] = 10

    fig = pplt.figure()
    ax = fig.add_subplot( 111, projection=ax3.Axes3D.name )

    xx = [ p.x for p in mesh.nodes_map.values() ]
    yy = [ p.y for p in mesh.nodes_map.values() ]
    zz = [ p.z for p in mesh.nodes_map.values() ]

    ax.set_xlabel('x')
    ax.set_ylabel('y')

    if add_labels:
        ax.plot( xx, yy, zz, 'o' )
        ax.legend( )

        cutoff_id = 1e9

        if len( mesh.nodes_map ) < 200:
            for p in mesh.nodes_map.values():
                if p.id < cutoff_id:
                    ax.text( p.x, p.y - 0.01, p.z, p.id, ha='right' )

        # plot the edges
        adj = mesh.get_adjacency_matrix()

        for entry in adj.items():
            edge = entry[0]

            if edge[0] < cutoff_id and edge[1] < cutoff_id:
                pt0 = mesh.nodes_map[ edge[0] ]
                pt1 = mesh.nodes_map[ edge[1] ]

                ax.plot( ( pt0.x, pt1.x ), ( pt0.y, pt1.y ), ( pt0.z, pt1.z ) )
    else:
        ax.plot_wireframe( xx, yy, zz )

    pplt.show()

    return


def run_plot_equatorial_points():

    import stomasimulator.stomata.stoma_config as sc
    cfg = sc.StomaConfig( cfg )

    cfg.print_detail()

    plot_equatorial_points( cfg )


def run_plot_mesh():
    import stomasimulator.stomata.stoma_config as sc

    cfg = default_cfg()
    # cfg = meckel_cfg()

    cfg.num_gcs = None

    cfg.print_detail()

    the_mesh = msh.build_stoma_mesh( cfg )

    plot_mesh( the_mesh )


def main():
    run_plot_equatorial_points()
    # run_plot_mesh()


if __name__ == '__main__':
    main()
