# ATHENA_Geant4

Geant4 model of the hadron endcap for the ATHENA EIC proposal. This is a small section of the endcap instead of the full version.

To use this model, clone and build this repository:
```
git clone https://github.com/rymilton/ATHENA_Geant4.git
cd ATHENA_Geant4
mkdir build && cd build
cmake ..
make install
````
If you make changes to the files, you'll have to `make install` again in the `build` directory. You can also install to another directory with
`cmake .. -DCMAKE_INSTALL_PREFIX=install` where `install` is some other directoy.

To run the simulation, go to your install directory and use `./ATHENA_Geometry -m mymac_WScFi.mac -t num_threads`,
where `num_threads` is the number of threads you want to use.
energy_loop.sh is also provided to generate multiple energies in one process.

For visualization, use `./ATHENA_Geometry`. Note to run simulations in this client, an output file name is still needed --
something like `/analysis/setFileName pi+_10GeV_5deg.root`.

During analysis, use the `eventID` variable to line up entries in different ntuples. An example analysis file is given with `Resolution.cpp`.
