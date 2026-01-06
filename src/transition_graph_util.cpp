// ====================================================================
// This file is part of PhaseTracer

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
// ====================================================================

#include <iostream>

#include "logger.hpp"
#include "transition_graph_util.hpp"

namespace TransitionGraph {

/**
 * Constructs a new phase from the input phase, with the position vector X being reflected about the axes defined in
 * reflectionIndices. For instance, if reflectionIndices = {0, 2} and phase.X[i] = (x0, x1, x2, x3), then
 * newPhase.X[i] = (-x0, x1, -x2, x3).
 *
 * @param phase - The phase explicitly stored in PhaseTracer, from which we derive the symmetric partners.
 * @param reflectionIndices - The axes in field space to reflect about for this phase. Note that this should be a
 *		subset of the overall symmetry of the potential. For instance, if the potential V(x0, x1) is symmetric in both
 *		field directions, then symmetryIndices should be {0}, {1}, and {0, 1} (in separate function calls) to handle
 *		all symmetric partners.
 * @param key - The key of the new phase.
 *
 * @return A new PhaseTracer::Phase identical to the input phase with the X property negated along the axes defined in
 *		reflectionIndices, and with the input key.
 */
PhaseTracer::Phase constructSymmetricPartnerPhase(const EffectivePotential::Potential &model, const PhaseTracer::Phase &phase, const std::vector<int> &reflectionIndices, int key) {
  PhaseTracer::Phase newPhase = phase;
  bool nonZero = false;

  // std::cout << std::endl;

  for (int i = 0; i < newPhase.X.size(); ++i) {
    for (int j = 0; j < reflectionIndices.size(); ++j) {
      // TODO: this criterion needs to scale with the scale of the fields.
      if (!nonZero && fabs(newPhase.X[i][reflectionIndices[j]]) > 0.01 * model.get_field_scale()) {
        nonZero = true;
      }

      newPhase.X[i][reflectionIndices[j]] *= -1;
    }
  }

  // A negative key signifies this is a redundant phase.
  newPhase.key = nonZero ? key : 0;

  return newPhase;
}

PhaseTracer::Transition constructSymmetricPartnerTransition(const PhaseTracer::Transition &transition, const PhaseTracer::Phase &false_phase, const PhaseTracer::Phase &true_phase, const Eigen::VectorXd &false_vacuum,
                                                            const Eigen::VectorXd &true_vacuum) {
  PhaseTracer::Transition newTransition = transition;

  // Can we avoid these copies?
  newTransition.false_phase = false_phase;
  newTransition.true_phase = true_phase;
  newTransition.false_vacuum = false_vacuum;
  newTransition.true_vacuum = true_vacuum;

  // Recalculate gamma since the distance between vacua has changed.
  // Unfortunately we can't call TransitionFinder's private method 'gamma'.
  newTransition.gamma = (true_vacuum - false_vacuum).norm() / transition.TC;

  return newTransition;
}

/**
 *	Generates all combinations of the symmetryIndices list, stored in a list of lists. O(n*N*(2^N)), where N is the
 * number of discrete symmetries and n is the number of fields.
 */
std::vector<std::vector<int>> generateReflectionLists(const std::vector<std::vector<int>> &symmetryIndices) {
  int N = symmetryIndices.size();
  int numSubLists = 1 << N;
  std::vector<std::vector<int>> reflectionLists;

  // This outer loop spans all combinations of the N symmetries.
  for (int i = 0; i < numSubLists; ++i) {
    std::vector<int> reflectionList = {};

    // For each combination, iterate over the N symmetries.
    for (int j = 0; j < N; ++j) {
      // Check if this combination should include this symmetry.
      if ((i >> j) & 1) {
        // Add all axes involved in this symmetry.
        for (int k = 0; k < symmetryIndices[j].size(); ++k) {
          reflectionList.push_back(symmetryIndices[j][k]);
        }
      }
    }

    reflectionLists.push_back(reflectionList);
  }

  return reflectionLists;
}

/** Returns the symmetric vacuum using the axes of reflection defined through reflectionIndices. */
Eigen::VectorXd getReflectedVacuum(const EffectivePotential::Potential &model, const Eigen::VectorXd &vacuum,
                                   const std::vector<int> &reflectionIndices, bool checkRedundancy) {
  Eigen::VectorXd newVacuum = vacuum;
  bool nonZero = false;

  if (checkRedundancy) {
    for (int i = 0; i < reflectionIndices.size(); ++i) {
      // TODO: this criterion needs to scale with the scale of the fields.
      if (!nonZero && fabs(newVacuum[reflectionIndices[i]]) > 0.001 * model.get_field_scale()) {
        nonZero = true;
      }

      newVacuum[reflectionIndices[i]] *= -1;
    }

    if (!nonZero) {
      newVacuum[0] = std::numeric_limits<double>::max();
    }
  } else {
    for (int i = 0; i < reflectionIndices.size(); ++i) {
      newVacuum[reflectionIndices[i]] *= -1;
    }
  }

  return newVacuum;
}

/** The comparison function for sorting phases by their keys. */
bool comparePhases(const PhaseTracer::Phase &a, const PhaseTracer::Phase &b) {
  return a.key < b.key;
}

/** The comparison function for sorting transitions by their false vacuum keys. */
bool compareTransitions(const PhaseTracer::Transition &a, const PhaseTracer::Transition &b) {
  return a.false_phase.key < b.false_phase.key;
}

// TODO: fix excessive copying?
void extractExplicitSymmetricPhasesAndTransitions(
    const EffectivePotential::Potential &model,
    const std::vector<PhaseTracer::Phase> &phases,
    const std::vector<PhaseTracer::Transition> &transitions,
    const std::vector<std::vector<int>> &symmetryIndices,
    std::vector<PhaseTracer::Phase> &out_symmetrisedPhases,
    std::vector<PhaseTracer::Transition> &out_symmetrisedTransitions) {
  // Generate all combinations of symmetryIndices.
  std::vector<std::vector<int>> reflectionLists = generateReflectionLists(symmetryIndices);

  // We will use this filtered list henceforth.
  std::vector<PhaseTracer::Transition> filteredTransitions;

  // Remove any symmetric transitions, or rather keep only the transitions between the explictly stored phases.
  for (int i = 0; i < transitions.size(); ++i) {
    if (transitions[i].key == 0) {
      filteredTransitions.push_back(transitions[i]);
    }
  }

  // Sort transitions by their false vacuum keys.
  std::sort(filteredTransitions.begin(), filteredTransitions.end(), compareTransitions);

  /*std::cout << "-------------------------------" << std::endl;
  std::cout << "Reflection lists:" << std::endl;
  for(int i = 0; i < reflectionLists.size(); ++i)
  {
          for(auto const& c : reflectionLists[i])
          {
                  std::cout << c << ' ';
          }

          std::cout << std::endl;
  }
  std::cout << "-------------------------------" << std::endl;*/

  // Prepare the symmetry lists for phases and transitions. There is a list for each phase and each transition. Each
  // list will be populated with the original phase or transition, then all symmetric partners of it.
  int numPhases = phases.size();
  int numTransitions = filteredTransitions.size();
  std::vector<std::vector<PhaseTracer::Phase>> phaseSymmetryList;
  std::vector<std::vector<PhaseTracer::Transition>> transitionSymmetryList;

  int maxKey = 0;

  for (int i = 0; i < numPhases; ++i) {
    std::vector<PhaseTracer::Phase> symmetryList;
    phaseSymmetryList.push_back(symmetryList);

    if (phases[i].key > maxKey) {
      maxKey = phases[i].key;
    }
  }

  for (int i = 0; i < numTransitions; ++i) {
    std::vector<PhaseTracer::Transition> symmetryList;
    transitionSymmetryList.push_back(symmetryList);
  }

  /*std::vector<int> keyToIndexMap(maxKey+1, 0);

  for(int i = 0; i < numPhases; ++i)
  {
          keyToIndexMap[phases[i].key] = i;
  }*/

  int key = phases.size();

  int transitionIndex = 0;

  // AHA! The issue is that using phase.key as an index is dangerous because if a phase is deleted (e.g. during
  // merging) then the largest key will exceed the bounds of the phases array.

  for (int i = 0; i < numPhases; ++i) {
    // phaseSymmetryList[phases[i].key].push_back(phases[i]);
    phaseSymmetryList[i].push_back(phases[i]);
    // Also add symmetric partners. Neglect the first 'empty' reflection.
    for (int j = 1; j < reflectionLists.size(); ++j) {
      PhaseTracer::Phase symmetricPhase = constructSymmetricPartnerPhase(model, phases[i], reflectionLists[j],
                                                                         key);

      // A key of zero signifies this reflection is redundant as the relevant field values are zero.
      // NOTE: the key is set specifically in constructSymmetricPartnerPhase to flag this.
      if (symmetricPhase.key > 0) {
        // phaseSymmetryList[phases[i].key].push_back(symmetricPhase);
        phaseSymmetryList[i].push_back(symmetricPhase);
        ++key;
      }
    }
  }

  for (int i = 0; i < numTransitions; ++i) {
    const PhaseTracer::Transition &transition = filteredTransitions[i];
    const PhaseTracer::Phase &falsePhase = transition.false_phase;
    const PhaseTracer::Phase &truePhase = transition.true_phase;

    // for(int j = 0; j < phaseSymmetryList[keyToIndexMap[falsePhase.key]].size(); ++j)
    for (int j = 0; j < phaseSymmetryList[falsePhase.key].size(); ++j) {
      const Eigen::VectorXd &symmetricFalseVacuum = getReflectedVacuum(model, transition.false_vacuum,
                                                                       reflectionLists[j], j > 0);

      // This flags that the reflected vacuum is redundant as the relevant field values are zero.
      if (symmetricFalseVacuum[0] == std::numeric_limits<double>::max()) {
        continue;
      }

      // for(int k = 0; k < phaseSymmetryList[keyToIndexMap[truePhase.key]].size(); ++k)
      for (int k = 0; k < phaseSymmetryList[truePhase.key].size(); ++k) {
        const Eigen::VectorXd &symmetricTrueVacuum = getReflectedVacuum(model, transition.true_vacuum,
                                                                        reflectionLists[k], k > 0);

        // This flags that the reflected vacuum is redundant as the relevant field values are zero.
        if (symmetricTrueVacuum[0] == std::numeric_limits<double>::max()) {
          continue;
        }

        /*transitionSymmetryList[i].push_back(constructSymmetricPartnerTransition(transition,
                phaseSymmetryList[keyToIndexMap[falsePhase.key]][j],
                phaseSymmetryList[keyToIndexMap[truePhase.key]][k], symmetricFalseVacuum, symmetricTrueVacuum));*/
        transitionSymmetryList[i].push_back(constructSymmetricPartnerTransition(transition,
                                                                                phaseSymmetryList[falsePhase.key][j], phaseSymmetryList[truePhase.key][k], symmetricFalseVacuum,
                                                                                symmetricTrueVacuum));
      }
    }
  }

  // std::cout << "================= SYMMETRIC PHASES ==================" << std::endl;

  for (int i = 0; i < numPhases; ++i) {
    for (int j = 0; j < phaseSymmetryList[i].size(); ++j) {
      out_symmetrisedPhases.push_back(phaseSymmetryList[i][j]);

      // std::cout << out_symmetrisedPhases[out_symmetrisedPhases.size()-1] << std::endl;
    }
  }

  // std::cout << "=============== SYMMETRIC TRANSITIONS ===============" << std::endl;

  for (int i = 0; i < numTransitions; ++i) {
    for (int j = 0; j < transitionSymmetryList[i].size(); ++j) {
      out_symmetrisedTransitions.push_back(transitionSymmetryList[i][j]);

      // std::cout << out_symmetrisedTransitions[out_symmetrisedTransitions.size()-1] << std::endl;
    }
  }

  // Sort so that the path finding algorithm correctly assumed phase index === phase key.
  // TODO: We can probably just relax that assumption in path finding...
  std::sort(out_symmetrisedPhases.begin(), out_symmetrisedPhases.end(), comparePhases);
}

PhaseStructureData extractPhaseStructureData(const EffectivePotential::Potential &model,
                                             const std::vector<PhaseTracer::Phase> &phases, const std::vector<PhaseTracer::Transition> &transitions,
                                             const std::vector<Eigen::VectorXd> &expectedLowTPhases, double Tmax, bool knownHighTPhase) {
  if (phases.size() == 0) {
    // std::cout << "No phases identified." << std::endl;
    return PhaseStructureData({-1}, {-1}, false, {});
  }

  std::vector<bool> isLowTemperaturePhase(phases.size(), false);

  /*if(transitions.size() == 0)
  {
          std::cout << "No transitions identified." << std::endl;
          return PhaseStructureData({-1}, {-1}, true, isLowTemperaturePhase);
  }*/

  // int EWVEVIndex = -1;
  // int globalMinimumIndex = -1;

  double minPotentialValue = std::numeric_limits<double>::max();
  std::vector<int> lowTPhaseIndices = {};

  double maxT = std::numeric_limits<double>::min();
  std::vector<int> highTPhaseIndices = {};

  for (int i = 0; i < phases.size(); ++i) {
    // At the start of this iteration, try to find if this phase is the high-T phase.
    if (phases[i].T.back() == Tmax) {
      // Only store one phase if we know the high-T phase.
      if (!knownHighTPhase || highTPhaseIndices.size() == 0) {
        highTPhaseIndices.push_back(i);
      } else if (phases[i].V.back() < phases[highTPhaseIndices[0]].V.back()) {
        highTPhaseIndices[0] = i;
      }
    }

    // For the rest of this iteration we check the phases at T=0.

    // Only check phases that exist at T=0. T.front() is the minimum temperature for the phase.
    if (phases[i].T.front() > 0) {
      continue;
    }

    /*Eigen::VectorXd zeroTPhase = phases[i].X.front();

    // Only worry about h > 0 phases, as the potential is symmetric in h and the EWVEV has h > 0.
    if(zeroTPhase[0] < -1.)
    {
            continue;
    }*/

    double zeroTPotential = phases[i].V.front();

    // Check if this is the deepest minimum (at T=0) found so far.
    if (zeroTPotential < minPotentialValue) {
      minPotentialValue = zeroTPotential;
      // globalMinimumIndex = i;

      // Clear the list of low-T phases since we found a lower energy phase, and add this phase to the list.
      lowTPhaseIndices = {i};

      // If we have found the EWVEV already, then we know it's not the global minimum at T=0.
      /*if(EWVEVIndex != -1)
      {
              std::cout << "A deeper minimum than the EWVEV was found." << std::endl;
              return PhaseStructureData(EWVEVIndex, {-1}, false);
      }*/
    }

    // TODO: this criterion needs to scale with the scale of the fields.
    /*if((zeroTPhase - EWVEV).squaredNorm() < 0.01)
    {
            // If a deeper minimum has already been found and it's not another phase coinciding with the EWVEV, then
            // we know the EWVEV is not the global minimum at T=0.
            if(globalMinimumIndex < i && globalMinimumIndex != EWVEVIndex)
            {
                    std::cout << "A deeper minimum than the EWVEV was found." << std::endl;
                    return PhaseStructureData(EWVEVIndex, {-1}, false);
            }

            EWVEVIndex = i;
    }*/
  }

  /*if(EWVEVIndex == -1)
  {
          std::cout << "No EWVEV at T=0." << std::endl;
          return PhaseStructureData(-1, highTPhaseIndices, false);
  }

  if(EWVEVIndex != globalMinimumIndex)
  {
          std::cout << "EWVEVIndex " << EWVEVIndex << " !=  globalMinimumIndex " << globalMinimumIndex;
          assert(EWVEVIndex == globalMinimumIndex && "Logic to detect global minima not equal to the EWVEV is flawed!");
  }*/

  for (int i = 0; i < lowTPhaseIndices.size(); ++i) {
    isLowTemperaturePhase[lowTPhaseIndices[i]] = true;
  }

  if (expectedLowTPhases.size() == 0) {
    return PhaseStructureData(lowTPhaseIndices, highTPhaseIndices, true, isLowTemperaturePhase);
  }

  bool valid = true;
  bool phaseMatchesExpected = true;

  // Check that each low-T phase was expected (i.e. it can be found in expectedLowTPhases). If this is not the case
  // for any of the phases, then the potential is labelled as invalid, since there is a low-T global minimum that was
  // not expected.
  for (int i = 0; i < lowTPhaseIndices.size(); ++i) {
    Eigen::VectorXd fieldVals = phases[lowTPhaseIndices[i]].X.front();

    // std::cout << "Checking low-T phase: (" << fieldVals[0] << ", " << fieldVals[1] << ")" << std::endl;

    for (int j = 0; j < expectedLowTPhases.size(); ++j) {
      phaseMatchesExpected = true;

      for (int k = 0; k < fieldVals.size(); ++k) {
        if (abs(fieldVals[k] - expectedLowTPhases[j][k]) > 0.01 * model.get_field_scale()) {
          phaseMatchesExpected = false;
          // std::cout << "fieldVals mismatch: k= " << k << "; (" << fieldVals[0] << ", " << fieldVals[1] <<
          //	") vs. (" << expectedLowTPhases[j][0] << ", " << expectedLowTPhases[j][1] << ")" << std::endl;
          break;
        }
      }

      if (phaseMatchesExpected) {
        // std::cout << "Phase matches expected: (" << fieldVals[0] << ", " << fieldVals[1] <<
        //		") vs. (" << expectedLowTPhases[j][0] << ", " << expectedLowTPhases[j][1] << ")" << std::endl;
        break;
      }
    }

    if (!phaseMatchesExpected) {
      LOG(debug) << "phase " << lowTPhaseIndices[i] << " does not match any expected low temperature phase." << std::endl;
      valid = false;
      break;
    }
  }

  // return PhaseStructureData(EWVEVIndex, highTPhaseIndices, true);
  return PhaseStructureData(lowTPhaseIndices, highTPhaseIndices, valid, isLowTemperaturePhase);
}

std::vector<Path> getTransitionPathsFromHighTPhase(const std::vector<Vertex> &vertices, int highTPhaseIndex,
                                                   const PhaseStructureData &phaseStructureData) {
  // Construct the initial frontier from all edges coming from the high temperature phase.
  int initialFrontierSize = vertices[highTPhaseIndex].edges.size();

  if (initialFrontierSize == 0) {
    LOG(debug) << "Initial frontier is empty." << std::endl;
    return {};
  }

  std::vector<FrontierNode> frontier;
  std::vector<Path> paths;
  std::vector<Path> terminalSubpaths;

  for (int i = 0; i < initialFrontierSize; ++i) {
    paths.push_back(Path(highTPhaseIndex));
    frontier.push_back(FrontierNode(highTPhaseIndex, i, i));
  }

  // Use a breadth-first search to find all edges from the high temperature phase to the EW VEV, using a frontier
  // approach. When the current end point of some path has multiple valid edges from it, we need to create new paths
  // that branch from this point, but are identical up to this point.
  while (frontier.size() > 0) {
    std::vector<FrontierNode> newFrontier;

    /*std::cout << "Frontier: {" << frontier[0];

    for(int i = 1; i < frontier.size(); ++i)
    {
            std::cout << ", " << frontier[i];
    }

    std::cout << "}" << std::endl;*/

    for (int i = 0; i < frontier.size(); ++i) {
      // std::cout << "Checking: " << frontier[i] << std::endl;

      const Edge &edge = vertices[frontier[i].vertexIndex].edges[frontier[i].edgeIndex];
      int vertexIndex = edge.toPhase;
      // Path& path = paths[frontier[i].pathIndex];

      double transitionTemperature = std::min(edge.temperature,
                                              paths[frontier[i].pathIndex].getCurrentTemperature());
      bool subcritical = transitionTemperature < edge.temperature;

      paths[frontier[i].pathIndex].extend(edge, subcritical, transitionTemperature);

      // This has been removed since we do care about extending the path further if there is a valid transition
      // away from the EW VEV. Of course, we need to relax conditions in extractPhaseStructureData for this to
      // be considered anyway.
      /*if(vertexIndex == phaseStructureData.EWVEVIndex)
      {
              // We're done with this path, so don't bother trying to extend it further.
              continue;
      }*/

      bool needToClonePath = false;

      // Check if the current vertex is one of the low temperature phases. If it is, we should add the current
      // path to the list of paths, prior to moving on to growing the path with outgoing edges from this vertex.

      int pathIndex = frontier[i].pathIndex;

      // std::cout << "Checking " << vertices[vertexIndex].edges.size() << " edges from vertex: " << vertexIndex
      //	<< std::endl;

      for (int j = 0; j < vertices[vertexIndex].edges.size(); ++j) {
        // if(paths[frontier[i].pathIndex].canUndergoTransition(vertices[vertexIndex].edges[j]))
        if (paths[frontier[i].pathIndex].canUndergoTransition(vertices, vertexIndex, j)) {
          if (needToClonePath) {
            // Call copy constructor.
            Path newPath = paths[frontier[i].pathIndex];
            paths.push_back(newPath);

            // This is merely assignment, and will not call the copy constructor again.
            // path = newPath;
            pathIndex = int(paths.size()) - 1;
          } else {
            needToClonePath = true;
          }

          // newFrontier.push_back(FrontierNode(vertices[vertexIndex].edges[j], path));
          newFrontier.push_back(FrontierNode(vertexIndex, j, pathIndex));
        }
      }

      if (needToClonePath) {
        /*for(int j = 0; j < phaseStructureData.lowTPhaseIndices.size(); ++j)
        {
                if(phaseStructureData.lowTPhaseIndices[j] == frontier[i].vertexIndex)
                {
                        // Create a copy of the current path.
                        Path newPath = path;

                        // The path has already been extended past the current vertex, so remove the last edge.
                        newPath.phases.pop_back();
                        newPath.temperatures.pop_back();
                        newPath.transitions.pop_back();

                        terminalSubpaths.push_back(newPath);

                        break;
                }
        }*/

        if (phaseStructureData.isLowTemperaturePhase[frontier[i].vertexIndex]) {
          /*std::cout << "At a low temperature phase, checking for terminal subpath..." << std::endl;

          std::cout << "Current terminal subpaths:" << std::endl;

          for(int j = 0; j < terminalSubpaths.size(); ++j)
          {
                  std::cout << j << ": " << terminalSubpaths[j] << std::endl;
          }

          std::cout << "Current path is:" << paths[frontier[i].pathIndex] << std::endl;*/

          // Check if this terminal subpath has already been added.
          bool bDuplicate = false;

          for (int j = 0; j < terminalSubpaths.size(); ++j) {
            // NOTE: the current path has one extra edge (the recently extended edge) that we will pop off
            // to obtain the terminal subpath. In the following, when we refer to the current path we will
            // neglect this last element/edge.

            // If the paths aren't the same size or don't end at the same low temperature phase, they are
            // not duplicate paths.
            if (terminalSubpaths[j].phases.size() != paths[frontier[i].pathIndex].phases.size() - 1 || terminalSubpaths[j].phases.back() != paths[frontier[i].pathIndex].phases.end()[-2]) {
              continue;
            }

            // This terminal subpath has the right properties to be a duplicate of the current path. We have
            // to explicitly match each phase along the paths to check if they're identical.
            bool mismatch = false;

            // Check transitions rather than phases. If the potential has phases that oscillate in their
            // relative energies, phases cannot be used to uniquely identify a path, but transitions can.
            for (int k = 0; k < terminalSubpaths[j].transitions.size(); ++k) {
              if (terminalSubpaths[j].transitions[k].transitionIndex !=
                  paths[frontier[i].pathIndex].transitions[k].transitionIndex) {
                mismatch = true;
                break;
              }
            }

            if (!mismatch) {
              // std::cout << "Detected as a duplicate of terminal subpath: " << j << std::endl;
              bDuplicate = true;
              break;
            }
          }

          if (!bDuplicate) {
            // std::cout << "Not a duplicate!" << std::endl;

            // Create a copy of the current path.
            Path newPath = paths[frontier[i].pathIndex];

            // The path has already been extended past the current vertex, so remove the last edge.
            newPath.phases.pop_back();
            // newPath.temperatures.pop_back();
            newPath.transitions.pop_back();

            terminalSubpaths.push_back(newPath);
          }
        }
      }
    }

    // Creates a copy of newFrontier.
    frontier = newFrontier;
  }

  std::vector<Path> validPaths;

  // std::cout << "Extracting valid paths..." << std::endl;

  for (int i = 0; i < paths.size(); ++i) {
    /*std::cout << "Path<" << i << ">" << std::endl;
    std::cout << "End phase: " << paths[i].phases.back() << std::endl;
    std::cout << "Is valid: " << phaseStructureData.isLowTemperaturePhase[paths[i].phases.back()] << std::endl;*/

    if (phaseStructureData.isLowTemperaturePhase[paths[i].phases.back()]) {
      validPaths.push_back(paths[i]);
    }
  }

  if (terminalSubpaths.size() > 0) {
    validPaths.insert(validPaths.end(), terminalSubpaths.begin(), terminalSubpaths.end());
  }

  LOG(debug) << "Transition paths from high-T phase:" << std::endl;
  LOG(debug) << "=============================" << std::endl;
  for (int i = 0; i < validPaths.size(); ++i) {
    LOG(debug) << validPaths[i] << std::endl;
  }
  LOG(debug) << "=============================" << std::endl;

  return validPaths;
}

std::vector<Path> getTransitionPaths(const std::vector<PhaseTracer::Phase> &phases,
                                     const std::vector<PhaseTracer::Transition> &transitions, const PhaseStructureData &phaseStructureData) {
  LOG(debug) << "TransitionGraph::getTransitionPaths" << std::endl;

  // std::vector<std::vector<Path>> paths(phaseStructureData.highTPhaseIndices.size());
  std::vector<Path> paths;

  // If there are no transitions, but one of the high temperature phases becomes a global minimum at T=0, then
  // we have a valid 'path' between the phases.
  for (int i = 0; i < phaseStructureData.highTPhaseIndices.size(); ++i) {
    for (int j = 0; j < phaseStructureData.lowTPhaseIndices.size(); ++j) {
      // if(phaseStructureData.highTPhaseIndices[i] == phaseStructureData.EWVEVIndex)
      if (phaseStructureData.highTPhaseIndices[i] == phaseStructureData.lowTPhaseIndices[j]) {
        // paths[i].push_back(Path(phaseStructureData.highTPhaseIndices[i]));
        paths.push_back(Path(phaseStructureData.highTPhaseIndices[i]));
        break;
      }
    }
  }

  if (transitions.size() == 0) {
    return paths;
  }

  // std::vector<Vertex> vertices(phases.size());
  std::vector<Vertex> vertices;

  // Construct the list of vertices from the phases.
  for (int i = 0; i < phases.size(); ++i) {
    // vertices[i] = Vertex(i);
    // vertices.push_back(Vertex(phases[i].key));
    vertices.push_back(Vertex(i));
  }

  // Construct the list of edges from the transitions.
  for (int i = 0; i < transitions.size(); ++i) {
    // std::cout << i << " " << transitions[i].id << std::endl;
    Edge edge(transitions[i].false_phase.key, transitions[i].true_phase.key, transitions[i].TC, transitions[i].id);
    vertices[edge.fromPhase].addEdge(edge);
  }

  /*std::cout << "============================================" << std::endl;
  std::cout << "Graph:" << std::endl;

  for(int i = 0; i < vertices.size(); ++i)
  {
          std::cout << vertices[i] << std::endl;
  }

  std::cout << "============================================" << std::endl;*/

  //::vector<Path> paths;

  for (int i = 0; i < phaseStructureData.highTPhaseIndices.size(); ++i) {
    LOG(debug) << "Finding transition paths from high-temperature phase: " << phaseStructureData.highTPhaseIndices[i] << std::endl;
    std::vector<Path> newPaths = getTransitionPathsFromHighTPhase(vertices, phaseStructureData.highTPhaseIndices[i],
                                                                  phaseStructureData);
    paths.insert(paths.end(), newPaths.begin(), newPaths.end());
  }

  return paths;
}

/**
 * knownHighTPhase=true means we know the global minimum at T=Tmax is the phase the Universe is in at T=Tmax.
 * knownHighTPhase=false means we don't know which of the phases at T=Tmax the Universe is in at T=Tmax.
 * This boolean variable is equivalent to whether we have sampled the potential at high enough temperatures.
 */
std::vector<Path> getPhaseHistory(const PhaseTracer::TransitionFinder &tf, bool knownHighTPhase) {
  std::vector<PhaseTracer::Phase> symmetrisedPhases;
  std::vector<PhaseTracer::Transition> symmetrisedTransitions;
  const EffectivePotential::Potential &model = tf.pf.get_potential(); 
  extractExplicitSymmetricPhasesAndTransitions(model, tf.pf.get_phases(), tf.get_transitions(),
                                               model.get_symmetry_axes(), symmetrisedPhases, symmetrisedTransitions);

  PhaseStructureData phaseStructureData = extractPhaseStructureData(model, symmetrisedPhases, symmetrisedTransitions,
                                                                    model.get_low_t_phases(),  tf.pf.get_t_high(), knownHighTPhase);

  if (!phaseStructureData.validAtZeroT) {
    LOG(debug) << "Phase structure is invalid at T=0. It does not describe our Universe.";
    return {};
  }

  //----------------------------------------------------------------------------------------------------
  // Report the global minima at T=Tmax and T=0.
  // std::cout << "Phase structure is valid at T=0." << std::endl;
  // std::cout << phaseStructureData << std::endl;
  /*std::cout << "Phase " << phaseStructureData.EWVEVIndex << " is the EW VEV phase at T=0." << std::endl;
  std::cout << "Phases ";

  for(int i = 0; i < phaseStructureData.highTPhaseIndices.size(); ++i)
  {
          std::cout << phaseStructureData.highTPhaseIndices[i] << " ";
  }

  if(phaseStructureData.highTPhaseIndices.size() == 1)
  {
          std::cout << "is the high temperature phase at T=";
  }
  else
  {
          std::cout << "are the high temperature phases at T=";
  }

  std::cout << pf.get_t_high() << "." << std::endl;*/

  //----------------------------------------------------------------------------------------------------
  // Report transition paths.
  // std::cout << std::endl;
  // std::cout << "Finding all transition paths from high temperature phases to low temperature phases..."
  //	<< std::endl;

  std::vector<Path> paths = getTransitionPaths(symmetrisedPhases, symmetrisedTransitions,
                                               phaseStructureData);

  return paths;

  /*int numPaths = 0;

  for(int i = 0; i < paths.size(); ++i)
  {
          for(int j = 0; j < paths[i].size(); ++j)
          {
                  numPaths += 1;
          }
  }

  std::cout << "Found " << numPaths << " paths:" << std::endl;

  int pathID = 1;
  //bool constrainedLowTPhases = model.get_low_t_phases().size() > 0;

  for(int i = 0; i < paths.size(); ++i)
  {
          for(int j = 0; j < paths[i].size(); ++j)
          {
                  //std::cout << (valid ? "[Valid]  " : "[Invalid]") << " Path " << pathID++  << ": " << paths[i][j]
                  //	<< std::endl;
                  std::cout << "Path " << pathID++ << ": " << paths[i][j] << std::endl;
                  //std::cout << paths[i] << std::endl;
          }
  }

  std::cout << std::endl;*/
}

} // namespace TransitionGraph
