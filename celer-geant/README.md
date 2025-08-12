Celeritas-Geant4 offloading app for validation
==============================================

Dependencies
------------

- Geant4 v11 or newer with `GEANT4_USE_GDML=ON` and
  `GEANT4_BUILD_MULTITHREADED=ON`
- Celeritas v0.6 or newer with `CELERITAS_USE_Geant4=ON`
- ROOT

Build and run
-------------

```sh
$ mkdir build
$ cd build
$ cmake ..
$ make
$ ./celer-geant input.json
```
