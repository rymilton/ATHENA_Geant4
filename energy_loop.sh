#!/bin/bash
set -e

filename="../mymac_WScFi.mac"
num_threads=14

particle="e-"

num_events=50000
energies=(0.5 1 2 5 10 20 30)

for (( k=1; k<7; k++ ))
do

	sed -i "s/\/analysis\/setFileName .*/\/analysis\/setFileName ${particle}_${energies[k]}GeV/" $filename
	sed -i "s/\/gps\/particle .*/\/gps\/particle ${particle}/" $filename
	sed -i "s/\/gps\/ene\/mono .*/\/gps\/ene\/mono ${energies[k]} GeV/" $filename
	sed -i "s/\/run\/beamOn .*/\/run\/beamOn ${num_events}/" $filename
	make
	echo "Working on energy ${energies[k]}"
	outfile="${particle}_gev${energies[k]}.log"
	echo "./ATHENA_Geometry -m mymac_WScFi.mac -t ${num_threads}"
	time ./ATHENA_Geometry -m mymac_WScFi.mac -t ${num_threads}
done
