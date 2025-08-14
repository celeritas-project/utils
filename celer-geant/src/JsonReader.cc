//------------------------------- -*- C++ -*- -------------------------------//
// Copyright Celeritas contributors: see top-level COPYRIGHT file for details
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file celer-geant/src/JsonReader.cc
//---------------------------------------------------------------------------//
#include "JsonReader.hh"

#include <fstream>
#include <corecel/Assert.hh>

//---------------------------------------------------------------------------//
//! Singleton declaration.
static JsonReader* json_reader_singleton{nullptr};

//---------------------------------------------------------------------------//
/*!
 * Construct singleton with JSON filename.
 */
void JsonReader::Construct(char const* json_filename)
{
    CELER_VALIDATE(!json_reader_singleton, << "JsonReader already constructed");
    json_reader_singleton = new JsonReader(json_filename);
}

//---------------------------------------------------------------------------//
/*!
 * Get JSON parser.
 *
 * \note \c JsonReader::construct(filename) must be called to construct it.
 */
nlohmann::json& JsonReader::Instance()
{
    CELER_VALIDATE(json_reader_singleton,
                   << "JsonReader not constructed. Initialize it by calling "
                      "JsonReader::construct(filename).");
    return json_reader_singleton->json_;
}

//---------------------------------------------------------------------------//
// PRIVATE
//---------------------------------------------------------------------------//
/*!
 * Construct from JSON filename.
 */
JsonReader::JsonReader(char const* json_filename)
{
    json_ = nlohmann::json::parse(std::ifstream(json_filename));
    CELER_VALIDATE(!json_.is_null(),
                   << "'" << json_filename << "' is not a valid input");
}
