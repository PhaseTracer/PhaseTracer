#pragma once

#include <vector>
#include <cmath>
#include <stdexcept>
#include <string>
#include <map>
#include <nlohmann/json.hpp>
#include <filesystem>
#include <random>
#include <algorithm>
#include "phasetracer.hpp"

using json = nlohmann::json;

namespace wallGoWrapper
{

/*
  Creates json object for wallGo to read in
*/
inline json
wallGoJSON(PhaseTracer::Transition& trans, const std::map<std::string, double>& param_map, const double& Tnuc) 
{

  // exact minima not needed, so use Tc values
  double h_min = abs(trans.true_vacuum[0]);
  double s_min = abs(trans.false_vacuum[1]);

  json input = {
    {"parameters", param_map},
    {"temperature", Tnuc},
    {"phase1", {0.0, s_min}},
    {"phase2", {h_min, 0.0}}
  };

  return input;
}

/*
  Write json object to file for wallGo to read in
*/
inline void
writeWallGoInputToJSON(const std::string& filename, const json& input) 
{
  std::ofstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error("Could not open file: " + filename);
  }
  file << input.dump(4);
  file.close();
}

/*
    Delete WallGoTempFile
*/
inline void
removeWallGoJSONFile(const std::string& filename)
{
    if(std::filesystem::exists(filename))
    {
        std::filesystem::remove(filename);
    }
}

inline std::string
getRandomFileName(const size_t length)
{
    const std::string hex_chars = "0123456789ABCDEF";
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distrib(0, hex_chars.length() - 1);

    std::string random_hex_string;
    random_hex_string.reserve(length + 5);

    for (size_t i = 0; i < length; ++i) {
        random_hex_string += hex_chars[distrib(gen)];
    }
    random_hex_string += ".json";

    return random_hex_string;
}

inline std::string 
execPythonScript(const std::string& cmd) 
{
  std::array<char, 256> buffer;
  std::string result;

  std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd.c_str(), "r"), pclose);
  if (!pipe) throw std::runtime_error("Failed to run Python script");

  while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
    result += buffer.data();
  }

  return result;
}

inline std::vector<double>
getWallVelocity(PhaseTracer::Transition& trans, const std::map<std::string, double>& param_map, const double& Tnuc, const std::string& model_name = "xSM")
{
  const auto filename = getRandomFileName(10);
  const json wallGoInputs = wallGoJSON(trans, param_map, Tnuc);
  writeWallGoInputToJSON(filename, wallGoInputs);

  std::string cmd = "python3 example/TestDeepPhase/helperIncludes/wallgo/" + model_name + ".py " + filename;

  const auto output = execPythonScript(cmd);

  std::vector<double> wallVelocities;
  std::istringstream iss(output);
  std::string token;
  
  // Parse comma-separated values
  while (std::getline(iss, token, ',')) {
    // Trim whitespace and newlines
    token.erase(std::remove_if(token.begin(), token.end(), ::isspace), token.end());
    if (!token.empty()) {
      wallVelocities.push_back(std::stod(token));
    }
  }

  removeWallGoJSONFile(filename);

  return wallVelocities;
}

} // namespace wallGoWrapper
