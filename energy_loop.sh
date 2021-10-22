#!/bin/bash
set -e

filename="../mymac_WScFi.mac"
num_threads=12

particle="pi+"

num_events=5000
energies=(1 2 5 10 20 30 40 50 60 70 80 90 100)

for (( k=8; k<13; k++ ))
do

	sed -i "s/\/analysis\/setFileName .*/\/analysis\/setFileName ${particle}_${energies[k]}GeV/" $filename
	sed -i "s/\/gps\/particle .*/\/gps\/particle ${particle}/" $filename
	sed -i "s/\/gps\/ene\/mono .*/\/gps\/ene\/mono ${energies[k]} GeV/" $filename
	sed -i "s/\/run\/beamOn .*/\/run\/beamOn ${num_events}/" $filename
	make
	echo "Working on energy ${energies[k]}"
	echo "./ATHENA_Geometry -m mymac_WScFi.mac -t ${num_threads}"
	time ./ATHENA_Geometry -m mymac_WScFi.mac -t ${num_threads}
done
