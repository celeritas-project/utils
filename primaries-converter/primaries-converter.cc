//------------------------------- -*- C++ -*- -------------------------------//
// Copyright Celeritas contributors: see top-level COPYRIGHT file for details
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file primaries-converter.cc
//---------------------------------------------------------------------------//
#include <iostream>
#include <celeritas/ext/RootFileManager.hh>
#include <corecel/cont/Range.hh>
#include <corecel/io/Logger.hh>

#include "src/EventReader.hh"
#include "src/JsonEventReader.hh"
#include "src/RootEventWriter.hh"

//---------------------------------------------------------------------------//
//! Loop over events of a given reader and write to ROOT.
template<class T>
void convert(T reader, RootEventWriter& writer, size_t max_events)
{
    size_t count = 0;
    for ([[maybe_unused]] auto i : celeritas::range(reader.num_events()))
    {
        if (max_events && count >= max_events)
            break;

        writer(reader());
        ++count;
    }

    CELER_LOG(info) << "Wrote " << count << " event(s) to ROOT file.";
}

//---------------------------------------------------------------------------//
/*!
 * Convert a HepMC3 or jsonl event record file to a ROOT file.
 */
int main(int argc, char* argv[])
{
    if (argc == 1)
    {
        std::cout << "Usage: " << argv[0]
                  << " input.[hepmc3/jsonl] [-n max_events] [-o output.root]"
                  << std::endl;
        return EXIT_FAILURE;
    }

    std::string input = argv[1];
    size_t max_events = 0;
    std::string root_filename;

    for (int i = 2; i < argc; ++i)
    {
        std::string arg = argv[i];
        if (arg == "-n")
        {
            if (i + 1 >= argc)
            {
                std::cout << "Error: missing value for " << arg << std::endl;
                return EXIT_FAILURE;
            }
            max_events = std::stoull(argv[++i]);
        }
        else if (arg == "-o")
        {
            if (i + 1 >= argc)
            {
                std::cout << "Error: missing value for " << arg << std::endl;
                return EXIT_FAILURE;
            }
            root_filename = argv[++i];
        }
        else
        {
            std::cout << "Error: unrecognized argument '" << arg << "'"
                      << std::endl;
            return EXIT_FAILURE;
        }
    }

    std::string extension = input.substr(input.find_last_of(".") + 1);

    if (extension != "hepmc3" && extension != "jsonl")
    {
        std::cout << "Error: input file must have .hepmc3 or .json extension"
                  << std::endl;
        return EXIT_FAILURE;
    }

    // Create ROOT file (default: input basename + ".root")
    if (root_filename.empty())
    {
        root_filename = input.substr(0, input.find_last_of("."));
        root_filename += ".root";
    }
    auto sp_rfm
        = std::make_shared<celeritas::RootFileManager>(root_filename.c_str());
    sp_rfm->make_tree("primaries", "primaries");
    RootEventWriter write_to_root(sp_rfm);

    // Write primaries to ROOT
    (extension == "hepmc3")
        ? convert(EventReader(input), write_to_root, max_events)
        : convert(JsonEventReader(input), write_to_root, max_events);

    return EXIT_SUCCESS;
}
