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

Adding new histograms
---------------------

The addition of any new histogram is straightforward:
- Expand JSON with new histogram information
- Add histogram to `SDHistograms` and initialize it in `SDHistograms::Initialize`
- Fill histogram in `SensitiveDetector::ProcessHits`
- Write it to disk during `RootIO::Finalize` using `RIO_HIST_WRITE`
