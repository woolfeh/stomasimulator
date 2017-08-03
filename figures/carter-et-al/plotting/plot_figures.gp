########################################################################################################################
# To run, start gnuplot in the same directory as this file
# Run...
# > call 'plot_figures.gp' 1
# where the 1 will produce eps files from the plots

reset

# 0 - do not output, 1 - output plots to eps files
#
output_plot=0

if ( ARGC > 0 ) {
  output_plot=1
}

########################################################################################################################
## Terminal id and plot placement variables

term_id = 0

sz_w  = 960; sz_h  = 540
scr_x = 200; scr_y = 100
dx    =  20; dy    = 20

########################################################################################################################
# Name the files

f_baseline    = '../simulations/spence1986-large-1um-with-epi-p,1gc,d_1.0,v_1.0,tiVW,c1_2.7579,c2_1.6122,c5_2320.35,lm_1.00.stats.txt'
f_vwt         = '../simulations/spence1986-large-1um-with-epi-p-rounded-triangle,1gc,d_1.0,v_1.0,tiVW,c1_2.7579,c2_1.6122,c5_2320.35,lm_1.00.stats.txt'
f_fixed_poles = '../simulations/spence1986-large-1um-with-epi-p,manually-edited-fixed-poles,1gc,d_1.0,v_1.0,tiVW,c1_2.7579,c2_1.6122,c5_2320.35,lm_1.00.stats.txt'

f_vwt_thicker = '../simulations/spence1986-large-1um-with-epi-p-rounded-thicker,1gc,d_1.0,v_1.1,tiVW,c1_2.7579,c2_1.6122,c5_2320.35,lm_1.00.stats.txt'
f_vwt_thinner = '../simulations/spence1986-large-1um-with-epi-p-rounded-thinner,1gc,d_1.0,v_0.9,tiVW,c1_2.7579,c2_1.6122,c5_2320.35,lm_1.00.stats.txt'

########################################################################################################################
# Pressure vs. aperture for baseline and VWT models

set term wxt term_id size sz_w,sz_h position scr_x + dx*term_id, scr_y + dy*term_id
term_id = term_id + 1

set xtics  out nomirror
set ytics  out mirror 0,2

set mxtics 2
set mytics 2

set border 1 + 2 + 8 # Left (2), bottom (1) and right (8)

set key top left Left reverse maxrows 3 width 1 height 0.5 samplen 3

set xlabel 'Turgor pressure (MPa)'
set xrange  [0:5]

set ylabel  'Aperture ({/Symbol m}m)'
set yrange [0:14]

set style rect fc lt -1 fs solid 0.1 noborder
set obj rect from 0, 0 to 0.5, graph 1

p \
  f_baseline    i 0 u 1:'pore-width' w l lw 3 lc 1 t'Baseline', \
  f_vwt         i 0 u 1:'pore-width' w l lw 3 lc 2 t'VWT'

if ( output_plot ) {
  call './to_eps.gp' './fig2B.eps' 'cols=1'
}

reset

########################################################################################################################
# Pressure vs. aperture for baseline with thicker and thinner

set term wxt term_id size sz_w,sz_h position scr_x + dx*term_id, scr_y + dy*term_id
term_id = term_id + 1

set xtics  out nomirror
set ytics  out mirror 0,2

set mxtics 2
set mytics 2

set border 1 + 2 + 8 # Left (2), bottom (1) and right (8)

set key top left Left reverse maxrows 3 width 1 height 0.5 samplen 3

set xlabel 'Turgor pressure (MPa)'
set xrange  [0:5]

set ylabel 'Aperture ({/Symbol m}m)'
set yrange [0:16]

set style rect fc lt -1 fs solid 0.1 noborder
set obj rect from 0, 0 to 0.5, graph 1

p \
  f_vwt            i 0 u 1:'pore-width' w l lw 3 lc 2      t'VWT', \
  f_vwt_thinner    i 0 u 1:'pore-width' w l lw 2 lc 2 dt 3 t'VWT: ventral wall 10% thinner', \
  f_vwt_thicker    i 0 u 1:'pore-width' w l lw 2 lc 2 dt 2 t'VWT: ventral wall 10% thicker'


if ( output_plot ) {
  call './to_eps.gp' './fig2C.eps' 'cols=1'
}

reset

########################################################################################################################
# Pressure vs. aperture for baseline and fixed pole models

set term wxt term_id size sz_w,sz_h position scr_x + dx*term_id, scr_y + dy*term_id
term_id = term_id + 1

set xtics  out nomirror
set ytics  out mirror 0,2

set mxtics 2
set mytics 2

set border 1 + 2 + 8 # Left (2), bottom (1) and right (8)

set key top left Left reverse maxrows 3 width 1 height 0.5 samplen 3

set xlabel 'Turgor pressure (MPa)'
set xrange  [0:5]

set ylabel  'Aperture ({/Symbol m}m)'
set yrange [0:16]

set style rect fc lt -1 fs solid 0.1 noborder
set obj rect from 0, 0 to 0.5, graph 1

p \
  f_baseline    i 0 u 1:'pore-width' w l lw 3 lc 1 t'Baseline', \
  f_fixed_poles i 0 u 1:'pore-width' w l lw 3 lc 3 t'Fixed poles'

if ( output_plot ) {
  call './to_eps.gp' './fig2D.eps' 'cols=1'
}

reset

########################################################################################################################
# Pressure vs. stoma length for baseline, VWT and fixed poles

set term wxt term_id size sz_w,sz_h position scr_x + dx*term_id, scr_y + dy*term_id
term_id = term_id + 1

set xtics  out nomirror
set ytics  out mirror 0,2

set mxtics 2
set mytics 2

set border 1 + 2 + 8 # Left (2), bottom (1) and right (8)

set key top left Left reverse maxrows 3 width 1 height 0.5 samplen 3

set xlabel 'Turgor pressure (MPa)'
set xrange  [0:5]

set ylabel 'Complex length ({/Symbol m}m)'
set yrange [44:52]

set style rect fc lt -1 fs solid 0.1 noborder
set obj rect from 0, 0 to 0.5, graph 1

p \
  f_baseline    i 0 u 1:'stoma-length' w l lw 3 lc 1 t'Baseline', \
  f_vwt         i 0 u 1:'stoma-length' w l lw 3 lc 2 t'VWT', \
  f_fixed_poles i 0 u 1:'stoma-length' w l lw 3 lc 3 t'Fixed poles'


if ( output_plot ) {
  call './to_eps.gp' './fig3A.eps' 'cols=1'
}

reset

########################################################################################################################



