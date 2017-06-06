########################################################################################################################
# Figure 1 - example model parameters for the Vicia faba data of Spence et al. (1986)
########################################################################################################################
# To run, start gnuplot in the same directory as this file
# Run...
# > call 'plot_fig1.gp'

reset

########################################################################################################################
## Terminal id and plot placement variables

term_id = 0

sz_w  = 960; sz_h  = 540
scr_x = 200; scr_y = 100
dx    =  20; dy    = 20

########################################################################################################################
# Get the files

# main iso
f_iso = '../simulations/vicia-example,2gc,d_1.0,v_1.0,MR,c1_9.0000,c2_9.0000.stats.txt'

# iso with epidermal pressure
f_iso_p = '../simulations/vicia-example-epi-p,1gc,d_1.0,v_1.0,MR,c1_9.0000,c2_9.0000.stats.txt'

# main aniso
f_aniso = '../simulations/vicia-example,2gc,d_1.0,v_1.0,tiMR,c1_9.0000,c2_9.0000,c5_1000.00,lm_1.00.stats.txt'

# aniso with epidermalpressure
f_aniso_p = '../simulations/vicia-example-epi-p,1gc,d_1.0,v_1.0,tiMR,c1_9.0000,c2_9.0000,c5_1000.00,lm_1.00.stats.txt'

########################################################################################################################
# Pressure vs. aperture and GC volume (iso)

set term wxt term_id size sz_w,sz_h position scr_x + dx*term_id, scr_y + dy*term_id
term_id = term_id + 1

set xtics  out nomirror
set ytics  out nomirror 0,1
set y2tics out nomirror
set mytics 2
set my2tics 2

set border 1 + 2 + 8 # Left (2), bottom (1) and right (8)

set key top center Left reverse maxrows 3 width 1 height 0.5 samplen 3

set xlabel 'Turgor pressure (MPa)'

set xrange  [0:5]
set x2range [0:5]

set ylabel  'Aperture ({/Symbol m}m)'
set y2label 'Guard cell volume ({/Symbol m}m^3)'

set yrange  [0:3]
set y2range [4000:12000]

p \
  f_iso     i 0 u 1:'pore-width' w l dt 2 lw 3 lc 1 t'Aperture'             axes x1y1, \
  ''	    i 0 u 1:'gc0-volume' w l dt 2 lw 3 lc 7 t'Guard cell volume'    axes x2y2

reset

########################################################################################################################
# Pressure vs. aperture and GC volume (aniso)

set term wxt term_id size sz_w,sz_h position scr_x + dx*term_id, scr_y + dy*term_id
term_id = term_id + 1

set xtics  out nomirror
set ytics  out nomirror 0,1
set y2tics out nomirror 4500,200
set mytics 2
set my2tics 2

set border 1 + 2 + 8 # Left (2), bottom (1) and right (8)

set key top center Left reverse maxrows 3 width 1 height 0.5 samplen 3

set xlabel 'Turgor pressure (MPa)'

set xrange  [0:5]
set x2range [0:5]

set ylabel  'Aperture ({/Symbol m}m)'
set y2label 'Guard cell volume ({/Symbol m}m^3)'

set yrange  [0:8]
set y2range [4500:6100]

p \
  f_aniso   i 0 u 1:'pore-width' w l lw 3 lc 1 t'Aperture'          axes x1y1, \
  ''        i 0 u 1:'gc0-volume' w l lw 3 lc 7 t'Guard cell volume' axes x2y2

reset

########################################################################################################################
# Pressure vs. aperture and GC volume (iso + epidermal pressure)

set term wxt term_id size sz_w,sz_h position scr_x + dx*term_id, scr_y + dy*term_id
term_id = term_id + 1

set xtics  out nomirror
set ytics  out nomirror 0,1
set y2tics out nomirror
set mytics 2
set my2tics 2

set border 1 + 2 + 8 # Left (2), bottom (1) and right (8)

set key top center Left reverse maxrows 3 width 1 height 0.5 samplen 3

set xlabel 'Turgor pressure (MPa)'

set xrange  [0:5]
set x2range [0:5]

set ylabel  'Aperture ({/Symbol m}m)'
set y2label 'Guard cell volume ({/Symbol m}m^3)'

set yrange  [0:3]
set y2range [4000:12000]

set style rect fc lt -1 fs solid 0.15 noborder
set obj rect from 0, 0 to 0.5, graph 1

p \
  f_iso_p   i 0 u 1:'pore-width' w l dt 2 lw 3 lc 1 t'Aperture'             axes x1y1, \
  ''        i 0 u 1:'gc0-volume' w l dt 2 lw 3 lc 7 t'Guard cell volume'    axes x2y2

reset

########################################################################################################################
# Pressure vs. aperture and GC volume (aniso + epidermal pressure)

set term wxt term_id size sz_w,sz_h position scr_x + dx*term_id, scr_y + dy*term_id
term_id = term_id + 1

set xtics  out nomirror
set ytics  out nomirror 0,1
set y2tics out nomirror 4500,200
set mytics 2
set my2tics 2

set border 1 + 2 + 8 # Left (2), bottom (1) and right (8)

set key top center Left reverse maxrows 3 width 1 height 0.5 samplen 3

set xlabel 'Turgor pressure (MPa)'

set xrange  [0:5]
set x2range [0:5]

set ylabel  'Aperture ({/Symbol m}m)'
set y2label 'Guard cell volume ({/Symbol m}m^3)'

set yrange  [0:8]
set y2range [4500:6100]

set style rect fc lt -1 fs solid 0.15 noborder
set obj rect from 0, 0 to 0.5, graph 1

p \
  f_aniso_p i 0 u 1:'pore-width' w l lw 3 lc 1 t'Aperture'          axes x1y1, \
  ''        i 0 u 1:'gc0-volume' w l lw 3 lc 7 t'Guard cell volume' axes x2y2

reset

########################################################################################################################
