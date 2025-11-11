#pragma once

#include <vector>
#include <cmath>
#include <stdexcept>
#include <string>
#include <map>
#include <fstream>
#include <cctype>
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
removeJSONFile(const std::string& filename)
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

struct wallGoResults
{
  // LTE limit
  double vwLTE;

  // from solveWall
  double vw_flag;
  double vw;
  double vJ;
  double vw_err;

  // from solveWallDetonation
  double vw_det_flag;
  double vw_det;
  double vJ_det;
  double vw_det_err;
};

inline wallGoResults
parseWallGoOutput(const std::string& wallGoFilename)
{
  std::string filename = wallGoFilename;
  filename.erase(std::remove_if(filename.begin(), filename.end(), [](unsigned char c) { return std::isspace(c); }), filename.end());


  std::ifstream file(filename);
  if (!file.is_open()) 
  {
    throw std::runtime_error("Could not open WallGo output file: " + filename);
  }

  json data;
  try 
  {
    file >> data;
  } catch (const json::parse_error& e) 
  {
    throw std::runtime_error("Failed to parse WallGo JSON output: " + std::string(e.what()));
  }
  file.close();

  wallGoResults results;
  
  results.vwLTE = data.value("vwLTE", -1.0);
  
  results.vw_flag = data.value("vw_flag", -1.0);
  results.vw = data.value("vw", -1.0);
  results.vJ = data.value("vJ", -1.0);
  results.vw_err = data.value("vw_err", -1.0);
  
  results.vw_det_flag = data.value("vwDet_flag", -1.0);
  results.vw_det = data.value("vwDet", -1.0);
  results.vJ_det = data.value("vJDet", -1.0);
  results.vw_det_err = data.value("vwDet_err", -1.0);

  return results;
}

inline wallGoResults
getWallVelocity(PhaseTracer::Transition& trans, const std::map<std::string, double>& param_map, const double& Tnuc, const std::string& model_name = "xSM")
{
  const auto filename = getRandomFileName(10);
  const json wallGoInputs = wallGoJSON(trans, param_map, Tnuc);
  writeWallGoInputToJSON(filename, wallGoInputs);

  std::string cmd = "python3 example/TestDeepPhase/helperIncludes/wallgo/" + model_name + ".py " + filename;

  const auto output = execPythonScript(cmd);

  const auto wallGoOutput = parseWallGoOutput(output);

  removeJSONFile(filename);
  removeJSONFile(output);

  return wallGoOutput;
}

} // namespace wallGoWrapper
