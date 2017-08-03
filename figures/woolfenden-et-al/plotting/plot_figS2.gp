########################################################################################################################
## Figure S2 - Thicker walls vs. stiffer matrix
########################################################################################################################

reset

########################################################################################################################
## Terminal id and plot placement variables

term_id = 0

sz_w  =  960; sz_h  = 540
scr_x = 1200; scr_y = 100
dx    =   20; dy    = 20

########################################################################################################################
# Name the files

# main iso
f_iso = '../simulations/vicia-example,2gc,d_1.0,v_1.0,MR,c1_9.0000,c2_9.0000.stats.txt'

# iso with different wall thicknesses
f_iso_d1v2 = '../simulations/vicia-example-thicker-ventral,2gc,d_1.0,v_2.0,MR,c1_9.0000,c2_9.0000.stats.txt'
f_iso_d2v2 = '../simulations/vicia-example-thicker,2gc,d_2.0,v_2.0,MR,c1_9.0000,c2_9.0000.stats.txt'

# iso with stiffer cell wall matrix

f_iso_g0_2x = '../simulations/vicia-example,MR,sens,G0x2.0.stats.txt'
f_iso_g0_4x = '../simulations/vicia-example,MR,sens,G0x4.0.stats.txt'

# main aniso
f_aniso = '../simulations/vicia-example,2gc,d_1.0,v_1.0,tiMR,c1_9.0000,c2_9.0000,c5_1000.00,lm_1.00.stats.txt'

# aniso with different wall thicknesses
f_aniso_d1v2 = '../simulations/vicia-example-thicker-ventral,1gc,d_1.0,v_2.0,tiMR,c1_9.0000,c2_9.0000,c5_1000.00,lm_1.00.stats.txt'
f_aniso_d2v2 = '../simulations/vicia-example-thicker,1gc,d_2.0,v_2.0,tiMR,c1_9.0000,c2_9.0000,c5_1000.00,lm_1.00.stats.txt'

# aniso with stiffer cell wall matrix

f_aniso_g0_2x = '../simulations/vicia-example,tiMR,sens,G0x2.0.stats.txt'
f_aniso_g0_4x = '../simulations/vicia-example,tiMR,sens,G0x4.0.stats.txt'

########################################################################################################################
# Pressure vs. relative change in aperture (for different cell wall thicknesses)

set term wxt term_id size sz_w,sz_h position scr_x + dx*term_id, scr_y + dy*term_id
term_id = term_id + 1

set xtics out nomirror
set ytics out
set border 1 + 2 + 8 # Left (2), bottom (1) and right (8)

set key top left Left reverse height 0.5 samplen 3

set xlabel 'Turgor pressure (MPa)'
set xrange [0:5]

set ylabel 'Aperture ({/Symbol m}m)'
set yrange [0:8]
set ytics 0,1
set mytics 2

p \
  f_aniso       i 0 u 1:'pore-width' w l lc 1      lw 3 t'With CMFs', \
  f_aniso_d1v2  i 0 u 1:'pore-width' w l lc 6      lw 3 t'With CMFs: thicker ventral wall', \
  f_aniso_d2v2  i 0 u 1:'pore-width' w l lc 4      lw 3 t'With CMFs: thicker walls', \
  f_iso         i 0 u 1:'pore-width' w l lc 1 dt 2 lw 3 t'Isotropic', \
  f_iso_d1v2    i 0 u 1:'pore-width' w l lc 6 dt 2 lw 3 t'Isotropic: thicker ventral wall', \
  f_iso_d2v2    i 0 u 1:'pore-width' w l lc 4 dt 2 lw 3 t'Isotropic: thicker walls'

reset

########################################################################################################################
# Pressure vs. relative change in aperture (for different G_0)

set term wxt term_id size sz_w,sz_h position scr_x + dx*term_id, scr_y + dy*term_id
term_id = term_id + 1

set xtics out nomirror
set ytics out
set border 1 + 2 + 8 # Left (2), bottom (1) and right (8)

set key top left Left reverse height 0.5 samplen 3

set xlabel 'Turgor pressure (MPa)'
set xrange [0:5]

set ylabel 'Aperture ({/Symbol m}m)'
set yrange [0:8]
set ytics 0,1
set mytics 2

p \
  f_aniso       i 0 u 1:'pore-width' w l lc 1      lw 3 t'With CMFs', \
  f_aniso_g0_2x i 0 u 1:'pore-width' w l lc 6      lw 3 t'With CMFs: 2x stiffer matrix', \
  f_aniso_g0_4x i 0 u 1:'pore-width' w l lc 4      lw 3 t'With CMFs: 4x stiffer matrix', \
  f_iso         i 0 u 1:'pore-width' w l lc 1 dt 2 lw 3 t'Isotropic', \
  f_iso_g0_2x   i 0 u 1:'pore-width' w l lc 6 dt 2 lw 3 t'Isotropic: 2x stiffer matrix', \
  f_iso_g0_4x   i 0 u 1:'pore-width' w l lc 4 dt 2 lw 3 t'Isotropic: 4x stiffer matrix'

reset

########################################################################################################################
