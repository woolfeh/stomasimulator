########################################################################################################################
## Figure 4 - Apertures for Arabidopsis opening experiments
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

cw_data = '../data/cw_data.txt'

# Col-0  cw1
# irx8   cw16
# pmr5   cw17
# pmr6   cw18

########################################################################################################################
# aperture histogram

set term wxt term_id size sz_w,sz_h position scr_x + dx*term_id, scr_y + dy*term_id
term_id = term_id + 1

set border 1 + 2 + 8 # Left (2), bottom (1) and right (8)

set xtics out scale 0 nomirror
set ytics out

set style fill solid 1.00 border 0
set style histogram errorbars gap 2 lw 1
set style data histogram

set boxwidth 0.25

set grid ytics

set xtics offset 1.9,0
set ytics 0,1
set mytics 2

set xtics ( "Col-0" 0.8, "{/Helvetica-Italic irx8}" 1.8, "{/Helvetica-Italic pmr5}" 2.8, "{/Helvetica-Italic pmr6}" 3.8 )

set xlabel "Genotype"
set ylabel "Aperture ({/Symbol m}m)"

set key top left Left reverse samplen 2

plot [0.5:4.5][0:3] \
    cw_data i 0 u ($0-0.15):2:3 with boxerrorbars   ti "Initial" lc 1, \
    ''      i 1 u ($0+0.15):2:3 with boxerrorbars   ti "Final"   lt 1 lc 4, \
    ''      i 0 u ($0-0.15):($2+$3+0.15):4 with labels notitle, \
    ''      i 1 u ($0+0.15):($2+$3+0.15):4 with labels notitle

reset

########################################################################################################################
