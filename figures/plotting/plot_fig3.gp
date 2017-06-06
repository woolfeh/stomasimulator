########################################################################################################################
## Figure 3 - Strain-stiffening cell wall matrix model applied to the Vicia faba data of Spence et al. (1986)
########################################################################################################################

reset

########################################################################################################################
## Terminal id and plot placement variables

term_id = 0

sz_w  =  960; sz_h  = 540
scr_x = 1200; scr_y = 100
dx    =   20; dy    = 20

# Borders: Bottom (1), Left (2) and right (8)
border_BLR = 1 + 2 + 8

########################################################################################################################
# Name the files

# optimised tiVW simulation
f_aniso_vw = '../simulations/vicia-optimised,2gc,d_1.0,v_1.0,tiVW,c1_2.7579,c2_1.6122,c5_2320.35,lm_1.00.stats.txt'

# get the previously computed tiMR results
f_aniso_mr = '../simulations/vicia-optimised,2gc,d_1.0,v_1.0,tiMR,c1_4.4546,c2_4.4779,c5_1169.65,lm_1.00.stats.txt'

# data from Franks et al. (2001), "Guard cell volume and pressure measured concurrently..."
franks_2001 = '../data/franks-2001-concurrent.txt'

########################################################################################################################
# Aperture vs. pressure

set term wxt term_id size sz_w,sz_h position scr_x + dx*term_id, scr_y + dy*term_id
term_id = term_id + 1

set xtics out nomirror
set ytics out
set border border_BLR

set key top left Left reverse

set xlabel 'Turgor pressure (MPa)'
set xrange [0:5]

set ylabel 'Aperture ({/Symbol m}m)'
set yrange [0:14]
set ytics 0,2

p \
  franks_2001   i 0 u 1:'half-aperture' w p pt 13 lc -1      t'Franks et al. (2001)', \
  f_aniso_vw    i 0 u 1:'pore-width'    w l lw 3  lc 1      t'Stiffening cell wall matrix', \
  f_aniso_mr    i 0 u 1:'pore-width'    w l lw 2  lc 1 dt 4 t'Original cell wall matrix'

reset

########################################################################################################################

set term wxt term_id size sz_w,sz_h position scr_x + dx*term_id, scr_y + dy*term_id
term_id = term_id + 1

set xtics out nomirror
set ytics out
set border border_BLR

set key top left Left reverse

set xlabel 'Turgor pressure (MPa)'
set xrange [0:5]

set ylabel 'Guard cell volume ({/Symbol m}m^3)'
set yrange [4500:7000]
set ytics 4500,500

p \
  franks_2001   i 0 u 1:'gc-volume'  w p pt 13 lc -1      t'Franks et al. (2001)', \
  f_aniso_vw    i 0 u 1:'gc0-volume' w l lw 3  lc 7      t'Stiffening cell wall matrix', \
  f_aniso_mr    i 0 u 1:'gc0-volume' w l lw 2  lc 7 dt 4 t'Original cell wall matrix'

reset

########################################################################################################################
