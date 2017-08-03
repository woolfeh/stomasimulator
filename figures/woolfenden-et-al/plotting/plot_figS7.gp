########################################################################################################################
## Figure S7 - Cell wall parameter sensitivity for tiVW model
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

f_aniso_tivw = '../simulations/vicia-optimised,2gc,d_1.0,v_1.0,tiVW,c1_2.7579,c2_1.6122,c5_2320.35,lm_1.00.stats.txt'

f_sens_g0 = system( 'ls ../simulations/vicia-optimised,tiVW,sens,G0x*stats.txt' )
f_sens_c5 = system( 'ls ../simulations/vicia-optimised,tiVW,sens,C5x*stats.txt' )

########################################################################################################################
# Vary G0 / Aperture vs. pressure

set term wxt term_id size sz_w,sz_h position scr_x + dx*term_id, scr_y + dy*term_id
term_id = term_id + 1

set xtics out nomirror
set ytics out
set border 1 + 2 + 8 # Left (2), bottom (1) and right (8)

set key top left Left reverse opaque height 0.5

set xlabel 'Turgor pressure (MPa)'

set ylabel 'Aperture ({/Symbol m}m)'
set yrange [0:22]
set ytics 0,2
set mytics 2

p \
  for [ f in f_sens_g0 ] f  i 0 u 1:'pore-width' w l lw 2 lc 1 dt 2 notitle, \
  f_aniso_tivw              i 0 u 1:'pore-width' w l lw 3 lc 1      notitle

reset

########################################################################################################################
# Vary G0 / GC width vs. pressure

set term wxt term_id size sz_w,sz_h position scr_x + dx*term_id, scr_y + dy*term_id
term_id = term_id + 1

set xtics out nomirror
set ytics out
set border 1 + 2 + 8 # Left (2), bottom (1) and right (8)

set key top left Left reverse opaque height 0.5

set xlabel 'Turgor pressure (MPa)'

set ylabel 'Guard cell width ({/Symbol m}m)'
set ytics format '%.1f'

p [][14.9:15.6] \
  for [ f in f_sens_g0 ] f  i 0 u 1:'gc0-width' w l lw 2 lc 4 dt 2 notitle, \
  f_aniso_tivw              i 0 u 1:'gc0-width' w l lw 3 lc 4      notitle

reset

########################################################################################################################
# Vary G0 / Volume vs. pressure

set term wxt term_id size sz_w,sz_h position scr_x + dx*term_id, scr_y + dy*term_id
term_id = term_id + 1

set xtics out nomirror
set ytics out
set border 1 + 2 + 8 # Left (2), bottom (1) and right (8)

set key top left Left reverse opaque height 0.5 width 1.5

set xlabel 'Turgor pressure (MPa)'

set ylabel 'Guard cell volume ({/Symbol m}m^3)'
set ytics 0,1000
set mytics 2

p [][4000:14000] \
  for [ f in f_sens_g0 ] f  i 0 u 1:'gc0-volume' w l lw 2 lc 7 dt 2 notitle, \
  f_aniso_tivw              i 0 u 1:'gc0-volume' w l lw 3 lc 7      notitle

reset

########################################################################################################################
# Vary C5 / Aperture vs. pressure

set term wxt term_id size sz_w,sz_h position scr_x + dx*term_id, scr_y + dy*term_id
term_id = term_id + 1

set xtics out nomirror
set ytics out
set border 1 + 2 + 8 # Left (2), bottom (1) and right (8)

set key top left Left reverse

set xlabel 'Turgor pressure (MPa)'

set ylabel 'Aperture ({/Symbol m}m)'
set yrange [0:22] #14]
set ytics 0,2
set mytics 2

p \
  for [ f in f_sens_c5 ] f  i 0 u 1:'pore-width' w l lw 2 lc 1 dt 2 notitle, \
  f_aniso_tivw              i 0 u 1:'pore-width' w l lw 3 lc 1      notitle

reset

########################################################################################################################
# Vary C5 / guard cell width vs. pressure

set term wxt term_id size sz_w,sz_h position scr_x + dx*term_id, scr_y + dy*term_id
term_id = term_id + 1

set xtics out nomirror
set ytics out
set border 1 + 2 + 8 # Left (2), bottom (1) and right (8)

set key top left Left reverse

set xlabel 'Turgor pressure (MPa)'

set ylabel 'Guard cell width ({/Symbol m}m)'
set mytics 2

p [][14:22] \
  for [ f in f_sens_c5 ] f  i 0 u 1:'gc0-width' w l lw 2 lc 4 dt 2 notitle, \
  f_aniso_tivw              i 0 u 1:'gc0-width' w l lw 3 lc 4      notitle

reset

########################################################################################################################
# Vary C5 / Volume vs. pressure

set term wxt term_id size sz_w,sz_h position scr_x + dx*term_id, scr_y + dy*term_id
term_id = term_id + 1

set xtics out nomirror
set ytics out
set border 1 + 2 + 8 # Left (2), bottom (1) and right (8)

set key top left Left reverse

set xlabel 'Turgor pressure (MPa)'

set ylabel 'Guard cell volume ({/Symbol m}m^3)'
set ytics 0,1000
set mytics 2

p [][4000:14000] \
  for [ f in f_sens_c5 ] f  i 0 u 1:'gc0-volume' w l lw 2 lc 7 dt 2 notitle, \
  f_aniso_tivw              i 0 u 1:'gc0-volume' w l lw 3 lc 7      notitle

reset

########################################################################################################################
