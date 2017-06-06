########################################################################################################################
## Figure S6 - Strain-stiffening cell wall matrix model applied to Vicia faba
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

########################################################################################################################
## Aperture and volume

set term wxt term_id size sz_w,sz_h position scr_x + dx*term_id, scr_y + dy*term_id
term_id = term_id + 1

set xtics  out nomirror
set ytics  out nomirror
set y2tics out nomirror
set mytics 2
set my2tics 2
set border 1 + 2 + 8 # Left (2), bottom (1) and right (8)

set key top left Left reverse

set xlabel 'Turgor pressure (MPa)'

set xrange  [0:5]
set x2range [0:5]

set ylabel  'Aperture ({/Symbol m}m)'
set y2label 'Guard cell volume ({/Symbol m}m^3)'

p \
	f_aniso_tivw i 0 u 1:'pore-width' w l lw 3 lc 1 t 'Aperture'           axes x1y1, \
    f_aniso_tivw i 0 u 1:'gc0-volume' w l lw 3 lc 7 t 'Guard cell volume'  axes x2y2

reset

########################################################################################################################
## GC surface area and stoma length

set term wxt term_id size sz_w,sz_h position scr_x + dx*term_id, scr_y + dy*term_id
term_id = term_id + 1

set xtics  out nomirror
set ytics  out nomirror
set y2tics out nomirror
set mytics 2
set my2tics 2
set border 1 + 2 + 8 # Left (2), bottom (1) and right (8)

set key top left Left reverse

set xlabel 'Turgor pressure (MPa)'

set xrange  [0:5]
set x2range [0:5]

set ylabel  'Stoma length ({/Symbol m}m)'
set y2label 'Guard cell surface area ({/Symbol m}m^2)'

set yrange  [45:53]
set y2range [1600:2400]

p \
	f_aniso_tivw i 0 u 1:'stoma-length'        w l lw 3 lc 2 t 'Stoma length'              axes x1y1, \
	f_aniso_tivw i 0 u 1:'gc0-surface-area'    w l lw 3 lc 4 t 'Guard cell surface area'   axes x2y2

reset

########################################################################################################################
