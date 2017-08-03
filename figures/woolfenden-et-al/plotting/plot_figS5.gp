########################################################################################################################
## Figure S5 - Strain-stiffening cell wall matrix (no fibres)
########################################################################################################################

reset

########################################################################################################################
## Terminal id and plot placement variables

term_id = 0

sz_w  = 960; sz_h  = 540
scr_x = 200; scr_y = 100
dx    =  20; dy    = 20

# Borders: Bottom (1), Left (2) and right (8)
border_BLR = 1 + 2 + 8

########################################################################################################################
# Name the files

f_iso_vw = '../simulations/vicia-optimised,2gc,d_1.0,v_1.0,VW,c1_2.7579,c2_1.6122.stats.txt'

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
set y2range [4000:16000]

p \
  f_iso_vw  i 0 u 1:'pore-width' w l dt 2 lw 3 lc 1 t'Aperture'             axes x1y1, \
  ''	    i 0 u 1:'gc0-volume' w l dt 2 lw 3 lc 7 t'Guard cell volume'    axes x2y2

reset

########################################################################################################################
