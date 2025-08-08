Celeritas-Geant4 offloading app for validation
==============================================

Dependencies
------------

- Geant4 v11 or newer
- Celeritas v0.5 or newer with ``CELERITAS_USE_Geant4=ON``
- ROOT

Build and run
-------------

```sh
$ mkdir build
$ cd build
$ cmake ..
$ make
$ export CELER_DISABLE_PARALLEL=1 # if Celeritas is built with MPI
$ ./celer-geant
```
