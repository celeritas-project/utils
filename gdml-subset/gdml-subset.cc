//----------------------------------*-C++-*----------------------------------//
// Copyright 2024 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file gdml-subset.cc
//---------------------------------------------------------------------------//
#include <cstdlib>
#include <string>
#include <CLI/CLI.hpp>
#include <G4GDMLParser.hh>
#include <G4PhysicalVolumeStore.hh>
#include <G4VPhysicalVolume.hh>
#include <G4Version.hh>
#include <corecel/Assert.hh>
#include <corecel/Version.hh>
#include <corecel/cont/Range.hh>
#include <corecel/io/Join.hh>
#include <corecel/io/Logger.hh>
#include <geocel/GeantGdmlLoader.hh>
#include <geocel/GeantGeoUtils.hh>
#include <geocel/ScopedGeantExceptionHandler.hh>
#include <geocel/ScopedGeantLogger.hh>

using namespace celeritas;
namespace
{
//---------------------------------------------------------------------------//
/*!
 * Print the usage of the app if possible, returning success.
 */
bool print_usage(CLI::App const& cli, std::ostream& os)
{
    if (auto base_formatter
        = std::dynamic_pointer_cast<CLI::Formatter>(cli.get_formatter()))
    {
        auto usage = base_formatter->make_usage(&cli, std::string{});
        if (!usage.empty() && usage.back() == '\n')
        {
            usage.pop_back();
        }
        os << usage;
        return true;
    }
    return false;
}

//---------------------------------------------------------------------------//
//! Construct a failure message for celeritas apps
std::string failure_message(CLI::App const* cli, const CLI::Error& e)
{
    std::ostringstream os;
    os << cli->get_name() << ": ";
    if (print_usage(*cli, os))
    {
        // Usage printed successfully; now write the error
        os << e.what();
    }
    else
    {
        // No usage available> write default error message
        os << CLI::FailureMessage::simple(cli, e);
    }

    return std::move(os).str();
}
struct Args
{
    std::string input_file;
    std::string volume_name;
    int depth{0};
    std::string output_file;
};

//---------------------------------------------------------------------------//
void delete_daughters_after(G4LogicalVolume* lv, int depth)
{
    if (depth == 0)
    {
        // Delete daughters
        lv->ClearDaughters();
        return;
    }
    --depth;

    for (auto const i : celeritas::range(lv->GetNoDaughters()))
    {
        delete_daughters_after(lv->GetDaughter(i)->GetLogicalVolume(), depth);
    }
}

G4VPhysicalVolume* find_volume(std::string const& vol_name)
{
    CELER_EXPECT(!vol_name.empty());

    auto& pvs = *G4PhysicalVolumeStore::GetInstance();
    auto new_world = std::find_if(
        pvs.begin(), pvs.end(), [&vol_name](G4VPhysicalVolume* pv) {
            return pv && pv->GetName() == vol_name;
        });

    CELER_VALIDATE(
        new_world != pvs.end(),
        << "failed to find volume '" << vol_name << "': available names are "
        << celeritas::join(
               pvs.begin(), pvs.end(), ", ", [](G4VPhysicalVolume* pv) {
                   return pv ? pv->GetName() : "<NULL>";
               }));
    return *new_world;
}

//---------------------------------------------------------------------------//
void run(Args const& args)
{
    // Read geometry *without* stripping pointers
    G4VPhysicalVolume* world = [&args] {
        using namespace celeritas;
        GeantGdmlLoader::Options opts;
        opts.pointers = GeantGdmlLoader::PointerTreatment::ignore;
        opts.detectors = false;
        return GeantGdmlLoader(opts)(args.input_file).world;
    }();

    // Find volume
    if (args.volume_name.empty())
    {
        CELER_LOG(info) << "Using original world volume";
    }
    else
    {
        world = find_volume(args.volume_name);
    }

    // Trim insides
    delete_daughters_after(world->GetLogicalVolume(), args.depth);

    // Write output
    G4GDMLParser parser;
    parser.SetEnergyCutsExport(false);
    parser.SetSDExport(false);
    parser.SetOverlapCheck(false);
#if G4VERSION_NUMBER >= 1070
    parser.SetOutputFileOverwrite(true);
#endif

    parser.Write(args.output_file, world, /* append_pointers = */ false);
}

}  // namespace

int main(int argc, char* argv[])
{
    Args args;

    static CLI::App app;
    app.failure_message(failure_message);
    app.set_version_flag("--version,-v", celeritas::version_string);
    app.description("Extract a subset of a GDML geometry file");
    app.add_option("--world",
                   args.volume_name,
                   "Physical volume name (empty string for world)");
    app.add_option("--depth", args.depth, "Depth to preserve (0 for all)")
        ->check(CLI::NonNegativeNumber);
    app.add_option("input", args.input_file, "Input GDML file")
        ->required()
        ->check(CLI::ExistingFile);
    app.add_option("output", args.output_file, "Output GDML file")->required();

    try
    {
        app.parse(argc, argv);
    }
    catch (CLI::ParseError const& e)
    {
        if (e.get_exit_code() != EXIT_SUCCESS)
        {
            world_logger()({app.get_name(), 0}, LogLevel::critical)
                << e.get_name() << ": " << e.what();
            print_usage(app, std::clog);
        }
        return app.exit(e);
    }

    try
    {
        ScopedGeantLogger scoped_log_;
        ScopedGeantExceptionHandler scoped_exceptions_;
        run(args);
    }
    catch (std::exception const& e)
    {
        auto msg = world_logger()({app.get_name(), 0}, LogLevel::critical);

        if (!dynamic_cast<RuntimeError const*>(&e))
        {
            // Not a Celeritas runtime error: print exception type
            msg << "Error: ";
        }
        msg << e.what();
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
