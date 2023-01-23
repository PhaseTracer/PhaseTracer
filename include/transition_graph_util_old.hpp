#ifndef TRANSITIONGRAPH_TRANSITION_GRAPH_UTIL_HPP_INCLUDED
#define TRANSITIONGRAPH_TRANSITION_GRAPH_UTIL_HPP_INCLUDED

#include "transition_finder.hpp"
#include "phase_finder.hpp"

#include <assert.h>
#include <stdlib.h>
#include <limits>

namespace PhaseTracer
{
	struct Transition;
}

namespace TransitionGraph
{

struct PhaseStructureData
{
	std::vector<int> lowTPhaseIndices;
	std::vector<int> highTPhaseIndices;
	bool validAtZeroT;

	PhaseStructureData(std::vector<int> lowTPhaseIndices, std::vector<int> highTPhaseIndices, bool validAtZeroT)
	{
		this->lowTPhaseIndices = lowTPhaseIndices;
		this->highTPhaseIndices = highTPhaseIndices;
		this->validAtZeroT = validAtZeroT;
	}

	friend std::ostream& operator << (std::ostream& o, const PhaseStructureData& p)
	{
		o << "Valid at T=0: " << p.validAtZeroT << ", low-T phases: {";

		if(p.lowTPhaseIndices.size() > 0)
		{
			o << p.lowTPhaseIndices[0];

			for(int i = 1; i < p.lowTPhaseIndices.size(); ++i)
			{
				o << " " << p.lowTPhaseIndices[i];
			}
		}

		o << "}, high-T phases: {";

		if(p.highTPhaseIndices.size() > 0)
		{
			o << p.highTPhaseIndices[0];

			for(int i = 1; i < p.highTPhaseIndices.size(); ++i)
			{
				o << " " << p.highTPhaseIndices[i];
			}
		}

		o << "}";
		return o;
	}
};

struct Edge
{
	int fromPhase;
	int toPhase;
	double temperature;
	int transition;

	Edge(int fromPhase, int toPhase, double temperature, int transition)
	{
		this->fromPhase = fromPhase;
		this->toPhase = toPhase;
		this->temperature = temperature;
		this->transition = transition;
	}

	friend std::ostream& operator << (std::ostream& o, const Edge& e)
	{	
		o << e.fromPhase << " --(" << e.transition << ", T=" << e.temperature << ")--> " << e.toPhase;
		return o;
	}
};

struct Path
{
	std::vector<int> phases;
	// This will be in descending order, with one less element than phases. temperatures[i] will be the critical
	// temperature between phases[i] and phases[i+1].
	std::vector<double> temperatures;
	// Similar to temperatures, this will be the transition indices in descending order of the corresponding
	// temperatures.
	std::vector<int> transitions;

	Path(int startPhase)
	{
		phases.push_back(startPhase);
	}

	// Copy constructor.
	/*Path(const Path& otherPath)
	{
		// These vectors will be copied to new vectors, which is good since we don't want changes to otherPath's
		// vectors to change this path's vectors.
		phases = otherPath.phases;
		temperatures = otherPath.temperatures;
	}*/

	void extend(const Edge& edge)
	{
		assert(edge.fromPhase == phases.back() && "Attempted to extend the path using an edge that doesn't begin at the end of this path!");

		phases.push_back(edge.toPhase);
		temperatures.push_back(edge.temperature);
		transitions.push_back(edge.transition);
	}

	bool canUndergoTransition(const Edge& edge)
	{
		return edge.temperature <= temperatures.back();
	}

	friend std::ostream& operator << (std::ostream& o, const Path& p)
	{
		if(p.phases.size() == 0)
		{
			o << "<empty path>";
			return o;
		}
		
		o << p.phases[0];

		for(int i = 1; i < p.phases.size(); ++i)
		{
			o << " --(" << p.transitions[i-1] << ", T=" << p.temperatures[i-1] << ")--> " << p.phases[i];
		}

		return o;
	}

	std::string getStringForFileOutput()
	{
		std::string transitionIndices = "";

		if(transitions.size() == 0)
		{
			return transitionIndices;
		}

		transitionIndices += std::to_string(transitions[0]);

		for(int i = 1; i < transitions.size(); ++i)
		{
			transitionIndices += " " + std::to_string(transitions[i]);
		}

		return transitionIndices;
	}
};

struct Vertex
{
	int phase;
	std::vector<Edge> edges;

	Vertex(int phase)
	{
		this->phase = phase;
		this->edges = {};
	}

	void addEdge(Edge& edge)
	{
		assert(edge.fromPhase == phase && "Attempted to add an edge to a vertex that is not the start point of the edge!");
		edges.push_back(edge);
	}

	friend std::ostream& operator << (std::ostream& o, const Vertex& v)
	{
		o << "Vertex<" << v.phase << "> {";

		if(v.edges.size() == 0)
		{
			o << "}";
			return o;
		}

		// Can't do size-1 if size=0 since it's stored in a size_t type var, which cannot store negative numbers.
		// That's why we have to handle the size = 0 case earlier.
		for(int i = 0; i < v.edges.size()-1; ++i)
		{
			o << v.edges[i] << ", ";
		}

		if(v.edges.size() > 0)
		{
			o << v.edges.back();
		}

		o << "}";

		return o;
	}
};

struct FrontierNode
{
	int vertexIndex;
	int edgeIndex;
	int pathIndex;

	FrontierNode(int vertexIndex, int edgeIndex, int pathIndex)
	{
		this->vertexIndex = vertexIndex;
		this->edgeIndex = edgeIndex;
		this->pathIndex = pathIndex;
	}

	friend std::ostream& operator << (std::ostream& o, const FrontierNode& fn)
	{	
		o << "{v: " << fn.vertexIndex << ", e: " << fn.edgeIndex << ", p: " << fn.pathIndex << "}";

		return o;
	}
};

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
PhaseTracer::Phase constructSymmetricPartnerPhase(const PhaseTracer::Phase& phase, const std::vector<int>&
	reflectionIndices, int key)
{
	PhaseTracer::Phase newPhase = phase;
	bool nonZero = false;

	std::cout << std::endl;

	for(int i = 0; i < newPhase.X.size(); ++i)
	{
		for(int j = 0; j < reflectionIndices.size(); ++j)
		{
			// TODO: this criterion needs to scale with the scale of the fields.
			if(!nonZero && fabs(newPhase.X[i][reflectionIndices[j]]) > 0.01)
			{
				nonZero = true;
			}
			
			newPhase.X[i][reflectionIndices[j]] *= -1;
		}
	}

	// A negative key signifies this is a redundant phase.
	newPhase.key = nonZero ? key : 0;

	return newPhase;
}

PhaseTracer::Transition constructSymmetricPartnerTransition(const PhaseTracer::Transition& transition, const 
	PhaseTracer::Phase& false_phase, const PhaseTracer::Phase& true_phase, const Eigen::VectorXd& false_vacuum,
	const Eigen::VectorXd& true_vacuum)
{
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
std::vector<std::vector<int>> generateReflectionLists(const std::vector<std::vector<int>>& symmetryIndices)
{
	int N = symmetryIndices.size();
	int numSubLists = 1 << N;
	std::vector<std::vector<int>> reflectionLists;

	// This outer loop spans all combinations of the N symmetries.
	for(int i = 0; i < numSubLists; ++i)
	{
		std::vector<int> reflectionList = {};

		// For each combination, iterate over the N symmetries.
		for(int j = 0; j < N; ++j)
		{
			// Check if this combination should include this symmetry.
			if((i >> j) & 1)
			{
				// Add all axes involved in this symmetry.
				for(int k = 0; k < symmetryIndices[j].size(); ++k)
				{
					reflectionList.push_back(symmetryIndices[j][k]);
				}
			}
		}

		reflectionLists.push_back(reflectionList);
	}

	return reflectionLists;
}

/** Returns the symmetric vacuum using the axes of reflection defined through reflectionIndices. */
Eigen::VectorXd getReflectedVacuum(const Eigen::VectorXd& vacuum, const std::vector<int>&
	reflectionIndices, bool checkRedundancy)
{
	Eigen::VectorXd newVacuum = vacuum;
	bool nonZero = false;

	if(checkRedundancy)
	{
		for(int i = 0; i < reflectionIndices.size(); ++i)
		{
			// TODO: this criterion needs to scale with the scale of the fields.
			if(!nonZero && fabs(newVacuum[reflectionIndices[i]]) > 0.01)
			{
				nonZero = true;
			}
		
			newVacuum[reflectionIndices[i]] *= -1;
		}

		if(!nonZero)
		{
			newVacuum[0] = std::numeric_limits<double>::max();
		}
	}
	else
	{
		for(int i = 0; i < reflectionIndices.size(); ++i)
		{
			newVacuum[reflectionIndices[i]] *= -1;
		}
	}

	return newVacuum;
}

/** The comparison function for sorting phases by their keys. */
bool comparePhases(const PhaseTracer::Phase& a, const PhaseTracer::Phase& b)
{
	return a.key < b.key;
}

/** The comparison function for sorting transitions by their false vacuum keys. */
bool compareTransitions(const PhaseTracer::Transition& a, const PhaseTracer::Transition& b)
{
	return a.false_phase.key < b.false_phase.key;
}

// TODO: fix excessive copying?
void extractExplicitSymmetricPhasesAndTransitions(
	const std::vector<PhaseTracer::Phase>& phases,
	const std::vector<PhaseTracer::Transition>& transitions,
	const std::vector<std::vector<int>>& symmetryIndices,
	std::vector<PhaseTracer::Phase>& out_symmetrisedPhases,
	std::vector<PhaseTracer::Transition>& out_symmetrisedTransitions)
{
	// Generate all combinations of symmetryIndices.
	std::vector<std::vector<int>> reflectionLists = generateReflectionLists(symmetryIndices);

	// We will use this filtered list henceforth.
	std::vector<PhaseTracer::Transition> filteredTransitions;

	// Remove any symmetric transitions, or rather keep only the transitions between the explictly stored phases.
	for(int i = 0; i < transitions.size(); ++i)
	{
		if(transitions[i].key == 0)
		{
			filteredTransitions.push_back(transitions[i]);
		}
	}

	// Sort transitions by their false vacuum keys.
	std::sort(filteredTransitions.begin(), filteredTransitions.end(), compareTransitions);

	std::cout << "-------------------------------" << std::endl;
	std::cout << "Reflection lists:" << std::endl;
	for(int i = 0; i < reflectionLists.size(); ++i)
	{
		for(auto const& c : reflectionLists[i])
		{
			std::cout << c << ' ';
		}

		std::cout << std::endl;
	}
	std::cout << "-------------------------------" << std::endl;

	// Prepare the symmetry lists for phases and transitions. There is a list for each phase and each transition. Each
	// list will be populated with the original phase or transition, then all symmetric partners of it.
	int numPhases = phases.size();
	int numTransitions = filteredTransitions.size();
	std::vector<std::vector<PhaseTracer::Phase>> phaseSymmetryList;
	std::vector<std::vector<PhaseTracer::Transition>> transitionSymmetryList;

	for(int i = 0; i < numPhases; ++i)
	{
		std::vector<PhaseTracer::Phase> symmetryList;
		phaseSymmetryList.push_back(symmetryList);
	}

	for(int i = 0; i < numTransitions; ++i)
	{
		std::vector<PhaseTracer::Transition> symmetryList;
		transitionSymmetryList.push_back(symmetryList);
	}

	int key = phases.size();

	int transitionIndex = 0;

	for(int i = 0; i < numPhases; ++i)
	{
		phaseSymmetryList[phases[i].key].push_back(phases[i]);
		// Also add symmetric partners. Neglect the first 'empty' reflection.
		for(int j = 1; j < reflectionLists.size(); ++j)
		{
			PhaseTracer::Phase symmetricPhase = constructSymmetricPartnerPhase(phases[i], reflectionLists[j], key);

			// A key of zero signifies this reflection is redundant as the relevant field values are zero.
			// NOTE: the key is set specifically in constructSymmetricPartnerPhase to flag this.
			if(symmetricPhase.key > 0)
			{
				phaseSymmetryList[phases[i].key].push_back(symmetricPhase);
				++key;
			}
		}
	}

	for(int i = 0; i < numTransitions; ++i)
	{
		const PhaseTracer::Transition& transition = filteredTransitions[i];
		const PhaseTracer::Phase& falsePhase = transition.false_phase;
		const PhaseTracer::Phase& truePhase = transition.true_phase;

		for(int j = 0; j < phaseSymmetryList[falsePhase.key].size(); ++j)
		{
			const Eigen::VectorXd& symmetricFalseVacuum = getReflectedVacuum(transition.false_vacuum,
				reflectionLists[j], j > 0);

			// This flags that the reflected vacuum is redundant as the relevant field values are zero.
			if(symmetricFalseVacuum[0] == std::numeric_limits<double>::max())
			{
				continue;
			}
			
			for(int k = 0; k < phaseSymmetryList[truePhase.key].size(); ++k)
			{
				const Eigen::VectorXd& symmetricTrueVacuum = getReflectedVacuum(transition.true_vacuum,
					reflectionLists[k], k > 0);

				// This flags that the reflected vacuum is redundant as the relevant field values are zero.
				if(symmetricTrueVacuum[0] == std::numeric_limits<double>::max())
				{
					continue;
				}

				transitionSymmetryList[i].push_back(constructSymmetricPartnerTransition(transition,
					phaseSymmetryList[falsePhase.key][j], phaseSymmetryList[truePhase.key][k],
					symmetricFalseVacuum, symmetricTrueVacuum));
			}
		}
	}

	std::cout << "================= SYMMETRIC PHASES ==================" << std::endl;

	for(int i = 0; i < numPhases; ++i)
	{
		for(int j = 0; j < phaseSymmetryList[i].size(); ++j)
		{
			out_symmetrisedPhases.push_back(phaseSymmetryList[i][j]);
			
			std::cout << out_symmetrisedPhases[out_symmetrisedPhases.size()-1] << std::endl;
		}
	}
	std::cout << "=============== SYMMETRIC TRANSITIONS ===============" << std::endl;

	for(int i = 0; i < numTransitions; ++i)
	{
		for(int j = 0; j < transitionSymmetryList[i].size(); ++j)
		{
			out_symmetrisedTransitions.push_back(transitionSymmetryList[i][j]);
			
			std::cout << out_symmetrisedTransitions[out_symmetrisedTransitions.size()-1] << std::endl;
		}
	}

	// Sort so that the path finding algorithm correctly assumed phase index === phase key.
	// TODO: We can probably just relax that assumption in path finding...
	std::sort(out_symmetrisedPhases.begin(), out_symmetrisedPhases.end(), comparePhases);
}

PhaseStructureData extractPhaseStructureData(const std::vector<PhaseTracer::Phase>& phases,
	const std::vector<PhaseTracer::Transition>& transitions, const std::vector<Eigen::VectorXd>& expectedLowTPhases,
	double Tmax, bool knownHighTPhase)
{
	if(phases.size() == 0)
	{
		std::cout << "No phases identified." << std::endl;
		return PhaseStructureData({-1}, {-1}, false);
	}

	if(transitions.size() == 0)
	{
		std::cout << "No transitions identified." << std::endl;
		return PhaseStructureData({-1}, {-1}, false);
	}

	//int EWVEVIndex = -1;
	//int globalMinimumIndex = -1;

	double minPotentialValue = std::numeric_limits<double>::max();
	std::vector<int> lowTPhaseIndices = {};

	double maxT = std::numeric_limits<double>::min();
	std::vector<int> highTPhaseIndices = {};

	for(int i = 0; i < phases.size(); ++i)
	{
		// At the start of this iteration, try to find if this phase is the high-T phase.
		if(phases[i].T.back() == Tmax)
		{
			// Only store one phase if we know the high-T phase.
			if(!knownHighTPhase || highTPhaseIndices.size() == 0)
			{
				highTPhaseIndices.push_back(i);
			}
			else if(phases[i].V.back() < phases[highTPhaseIndices[0]].V.back())
			{
				highTPhaseIndices[0] = i;
			}
		}
		
		// For the rest of this iteration we check the phases at T=0.

		// Only check phases that exist at T=0. T.front() is the minimum temperature for the phase.
		if(phases[i].T.front() > 0)
		{
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
		if(zeroTPotential < minPotentialValue)
		{
			minPotentialValue = zeroTPotential;
			//globalMinimumIndex = i;

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

	if(expectedLowTPhases.size() == 0)
	{
		return PhaseStructureData(lowTPhaseIndices, highTPhaseIndices, true);
	}

	bool valid = true;
	bool phaseMatchesExpected = true;

	// Check that each low-T phase was expected (i.e. it can be found in expectedLowTPhases). If this is not the case
	// for any of the phases, then the potential is labelled as invalid, since there is a low-T global minimum that was
	// not expected.
	for(int i = 0; i < lowTPhaseIndices.size(); ++i)
	{
		phaseMatchesExpected = true;

		for(int j = 0; j < expectedLowTPhases.size(); ++j)
		{
			Eigen::VectorXd fieldVals = phases[i].X.front();

			for(int k = 0; k < fieldVals.size(); ++k)
			{
				if(abs(fieldVals[k] - expectedLowTPhases[j][k]) > 0.01)
				{
					phaseMatchesExpected = false;
					break;
				}
			}

			if(phaseMatchesExpected)
			{
				break;
			}
		}

		if(!phaseMatchesExpected)
		{
			valid = false;
			break;
		}
	}
	
	//return PhaseStructureData(EWVEVIndex, highTPhaseIndices, true);
	return PhaseStructureData(lowTPhaseIndices, highTPhaseIndices, valid);
}

std::vector<Path> getTransitionPathsFromHighTPhase(const std::vector<Vertex>& vertices, int highTPhaseIndex)
{
	// Construct the initial frontier from all edges coming from the high temperature phase.
	int initialFrontierSize = vertices[highTPhaseIndex].edges.size();

	if(initialFrontierSize == 0)
	{
		return {};
	}

	std::vector<FrontierNode> frontier;
	std::vector<Path> paths;

	for(int i = 0; i < initialFrontierSize; ++i)
	{
		paths.push_back(Path(highTPhaseIndex));
		frontier.push_back(FrontierNode(highTPhaseIndex, i, i));
	}
	
	// Use a breadth-first search to find all edges from the high temperature phase to the EW VEV, using a frontier
	// approach. When some the current end point of some path has multiple valid edges from it, we need to create new
	// paths that branch from this point, but are identical up to this current point.
	while(frontier.size() > 0)
	{
		std::vector<FrontierNode> newFrontier;

		for(int i = 0; i < frontier.size(); ++i)
		{
			const Edge& edge = vertices[frontier[i].vertexIndex].edges[frontier[i].edgeIndex];
			int vertexIndex = edge.toPhase;
			Path& path = paths[frontier[i].pathIndex];

			path.extend(edge);

			// This has been removed since we do care about extending the path further if there is a valid transition
			// away from the EW VEV. Of course, we need to relax conditions in extractPhaseStructureData for this to 
			// be considered anyway.
			/*if(vertexIndex == phaseStructureData.EWVEVIndex)
			{
				// We're done with this path, so don't bother trying to extend it further.
				continue;
			}*/

			bool foundValidEdge = false;

			int pathIndex = frontier[i].pathIndex;

			for(int j = 0; j < vertices[vertexIndex].edges.size(); ++j)
			{
				if(path.canUndergoTransition(vertices[vertexIndex].edges[j]))
				{
					if(foundValidEdge)
					{
						// Call copy constructor.
						Path newPath = path;
						paths.push_back(newPath);
						
						// This is merely assignment, and will not call the copy constructor again.
						path = newPath;
						pathIndex = paths.size()-1;
					}

					//newFrontier.push_back(FrontierNode(vertices[vertexIndex].edges[j], path));
					newFrontier.push_back(FrontierNode(vertexIndex, j, pathIndex));
					foundValidEdge = true;
				}
			}
		}
		
		// Creates a copy of newFrontier.
		frontier = newFrontier;
	}

	return paths;
}

std::vector<std::vector<Path>> getTransitionPaths(const std::vector<PhaseTracer::Phase>& phases,
	const std::vector<PhaseTracer::Transition>& transitions, const PhaseStructureData& phaseStructureData)
{
	if(transitions.size() == 0)
	{
		std::vector<std::vector<Path>> paths(phaseStructureData.highTPhaseIndices.size());

		// If there are no transitions, but one of the high temperature phases becomes a global minimum at T=0, then
		// we have a valid 'path' between the phases.
		for(int i = 0; i < phaseStructureData.highTPhaseIndices.size(); ++i)
		{
			for(int j = 0; j < phaseStructureData.lowTPhaseIndices.size(); ++j)
			{
				//if(phaseStructureData.highTPhaseIndices[i] == phaseStructureData.EWVEVIndex)
				if(phaseStructureData.highTPhaseIndices[i] == phaseStructureData.lowTPhaseIndices[j])
				{
					paths[i].push_back(Path(phaseStructureData.highTPhaseIndices[i]));
					break;
				}
			}
		}
		
		return paths;
	}
	
	//std::vector<Vertex> vertices(phases.size());
	std::vector<Vertex> vertices;

	// Construct the list of vertices from the phases.
	for(int i = 0; i < phases.size(); ++i)
	{
		//vertices[i] = Vertex(i);
		//vertices.push_back(Vertex(phases[i].key));
		vertices.push_back(Vertex(i));
	}
	
	// Construct the list of edges from the transitions.
	for(int i = 0; i < transitions.size(); ++i)
	{
		std::cout << i << " " << transitions[i].id << std::endl;
		Edge edge(transitions[i].false_phase.key, transitions[i].true_phase.key, transitions[i].TC, transitions[i].id);
		vertices[edge.fromPhase].addEdge(edge);
	}

	std::cout << "============================================" << std::endl;
	std::cout << "Graph:" << std::endl;

	for(int i = 0; i < vertices.size(); ++i)
	{
		std::cout << vertices[i] << std::endl;
	}

	std::cout << "============================================" << std::endl;

	std::vector<std::vector<Path>> paths;

	for(int i = 0; i < phaseStructureData.highTPhaseIndices.size(); ++i)
	{
		paths.push_back(getTransitionPathsFromHighTPhase(vertices, phaseStructureData.highTPhaseIndices[i]));
	}

	return paths;
}

/**
  * knownHighTPhase=true means we know the global minimum at T=Tmax is the phase the Universe is in at T=Tmax.
  * knownHighTPhase=false means we don't know which of the phases at T=Tmax the Universe is in at T=Tmax.
  * This boolean variable is equivalent to whether we have sampled the potential at high enough temperatures.
  */
std::vector<std::vector<Path>> getPhaseHistory(const PhaseTracer::PhaseFinder& pf, const
	PhaseTracer::TransitionFinder& tf, const EffectivePotential::Potential& model, bool knownHighTPhase)
{
	std::vector<PhaseTracer::Phase> symmetrisedPhases;
	std::vector<PhaseTracer::Transition> symmetrisedTransitions;

	extractExplicitSymmetricPhasesAndTransitions(pf.get_phases(), tf.get_transitions(), model.get_symmetry_axes(),
		symmetrisedPhases, symmetrisedTransitions);

	PhaseStructureData phaseStructureData = extractPhaseStructureData(symmetrisedPhases, symmetrisedTransitions,
		model.get_low_t_phases(), pf.get_t_high(), knownHighTPhase);

	if(!phaseStructureData.validAtZeroT)
	{
		std::cout << "Phase structure is invalid at T=0. It does not describe our Universe." << std::endl;
		return {};
	}

	//----------------------------------------------------------------------------------------------------
	// Report the global minima at T=Tmax and T=0.
	std::cout << "Phase structure is valid at T=0." << std::endl;
	std::cout << phaseStructureData << std::endl;
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
	std::cout << std::endl;
	std::cout << "Finding all transition paths from high temperature phases to low temperature phases..."
		<< std::endl;

	std::vector<std::vector<Path>> paths = getTransitionPaths(symmetrisedPhases, symmetrisedTransitions,
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

#endif
