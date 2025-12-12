Primaries input converter
=========================

Convert a HepMC3 or jsonl file with primary particles exported by Celeritas to a
Celeritas' compatible ROOT file. This is useful in cases where it is not trivial
to dump primaries into a ROOT file directly, such as Athena. The use of ROOT
substantially reduces file size and load time.

# Dependencies

- ROOT
- HepMC3
- nlohmann/json
- Celeritas

# Build and run

```sh
$ mkdir build; cd build
$ cmake ..
$ make
$ ./convert input.[hepmc3/jsonl]
```
