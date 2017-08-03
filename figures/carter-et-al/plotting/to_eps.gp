# Output plot to an eps file
# Usage:
# > call 'to_eps.gp' 'filename' {keep_key}
# where {keep_key} is only required if the key is wanted (default is to omit)

# Output will be written to current directory

# arg0 is script name
# arg1 is filename
# arg2 is a flag to keep the key (valid options are 'show-key=no' and 'show-key=yes')
# arg3 is the number of columns that the plot will span (valid options are 'cols=1' and 'cols=2'
# argc is number of parameters (excluding arg0)
#
if ( ARGC == 0 ) {
  print "Incorrect usage: call 'to_eps.gp' 'filename' 'show-key={yes*|no}' 'cols={1*|2}"
  print "where the * indicates the default"
  exit
}

show_key = 1
if ( ARGC >= 2 ) {
  if ( "".ARG2 eq "show-key=no" || "".ARG3 eq "show-key=no" ) {
    show_key = 0
  }
}

cols = 1
if ( ARGC >= 2 ) {
  if ( "".ARG2 eq "cols=2" || "".ARG3 eq "cols=2" ) {
    cols = 2
  }
}


# Guideline widths:
# 1 column figure: 20 picas (3.4 in)
# 2 column figure: 40 picas (6.7 in)

filename=''.ARG1

#set term epscairo size 25.6cm, 16cm enhanced color font 'Helvetica,20' lw 4

if ( cols == 1 ) {
  set term epscairo size 3.4in, 2.125in enhanced color fontscale 1 font 'Helvetica,6' lw 2
} else {
  set term epscairo size 6.7in, 4.1875in enhanced color fontscale 1 font 'Helvetica,6' lw 2
}

set output filename

if ( show_key ) {
  set key font 'Helvetica,6'
} else {
  unset key
}

replot

set output
set term pop

print '--> Wrote plot to '.filename

if ( ARGC == 2 ) {
  set key default
}
