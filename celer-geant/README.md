Celeritas-Geant4 offloading app for validation
==============================================

# Dependencies

- Geant4 v11 or newer with `GEANT4_USE_GDML=ON` and
  `GEANT4_BUILD_MULTITHREADED=ON`
- Celeritas v0.6 or newer with `CELERITAS_USE_Geant4=ON`
- ROOT

# Build and run

```sh
$ mkdir build
$ cd build
$ cmake ..
$ make
$ ./celer-geant input.json
```

# Adding new histograms

The addition of any new histogram is straightforward:
- Expand JSON with new histogram information
- `HistogramStore.hh`: Add histogram to `SDHistograms` and initialize it in
  `SDHistograms::Initialize`
- `SensitiveDetector.cc`: Fill histogram in `SensitiveDetector::ProcessHits`
- `RootIO.cc`: Write histogram to disk during `RootIO::Finalize` using
  `RIO_HIST_WRITE`
