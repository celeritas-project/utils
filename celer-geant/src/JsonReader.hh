//------------------------------- -*- C++ -*- -------------------------------//
// Copyright Celeritas contributors: see top-level COPYRIGHT file for details
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file celer-geant/src/JsonReader.hh
//---------------------------------------------------------------------------//
#pragma once

#include <string>
#include <nlohmann/json.hpp>

//---------------------------------------------------------------------------//
/*!
 * Singleton \e nlohmann/json parser.
 * Use \c JsonReader::Construct("input.json") to construct the singleton, and
 * \c JsonReader::Instance() to access it from any class.
 */
class JsonReader
{
  public:
    //! Construct by creating singleton from the json filename
    static void Construct(char const* json_filename);

    //! Instance singleton with json parser
    static nlohmann::json& Instance();

    //! Validate JSON entry
    static void Validate(nlohmann::json const& j, std::string name);

    //! Validate JSON histogram entry
    static void
    ValidateHistogram(nlohmann::json const& j, std::string hist_name);

  private:
    // Json parser
    nlohmann::json json_;

    // Construct with filename
    JsonReader(char const* json_filename);
};
