#!/usr/bin/env bash

# make sure:
# - the FEBIO variable below is set to the location of the FEBio executable
# - the Python environment for 'stomasimulator' is active

#FEBIO=<path/to/febio>/FEBio2
FEBIO=${HOME}/Applications/febio-apps/febio/current/bin/FEBio2

echo "Running the simulations..."

for feb_file in *.feb
do
	${FEBIO} -silent -i ${feb_file} &
done

wait

echo "Processing the simulation files..."

for xplt_file in *.xplt
do
	stomasimulator xplt ${xplt_file}
done

echo "Finished"
