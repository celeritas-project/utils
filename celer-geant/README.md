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

- Expand JSON with new histogram information.
- `RootDataStore.hh`: Add histogram to `SensDetData` and initialize it in
  `SensDetData::Initialize`.
- Fill histogram (usually via `SensitiveDetector::ProcessHits`).
- `RootIO.cc`: Write histogram to disk during `RootIO::Finalize` using
  `RIO_HIST_WRITE`.

# I/O
Since `RootIO` is a thread-local singleton that owns a `RootDataStore` object,
which maps all sensitive detector data, it can be used as a starting point to
expand it to a more comprehensive I/O system that manages more complex data and
different output types. 
