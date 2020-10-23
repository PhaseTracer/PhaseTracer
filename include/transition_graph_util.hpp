#ifndef TRANSITIONGRAPH_TRANSITION_GRAPH_UTIL_HPP_INCLUDED
#define TRANSITIONGRAPH_TRANSITION_GRAPH_UTIL_HPP_INCLUDED

#include "transition_finder.hpp"
#include "phase_finder.hpp"

#include <assert.h>
#include <stdlib.h>
#include <limits>

namespace TransitionGraph
{

struct PhaseStructureData
{
	int EWVEVIndex;
	int highTPhaseIndex;
	bool validAtZeroT;

	/*PhaseStructureData()
	{
		assert(false && "Should never be calling PhaseStructureData's default constructor!");
	}*/

	PhaseStructureData(int EWVEVIndex, int highTPhaseIndex, bool validAtZeroT)
	{
		this->EWVEVIndex = EWVEVIndex;
		this->highTPhaseIndex = highTPhaseIndex;
		this->validAtZeroT = validAtZeroT;
	}
};

struct Edge
{
	int fromPhase;
	int toPhase;
	double temperature;

	/*Edge()
	{
		assert(false && "Should never be calling Edge's default constructor!");
	}*/

	Edge(int fromPhase, int toPhase, double temperature)
	{
		this->fromPhase = fromPhase;
		this->toPhase = toPhase;
		this->temperature = temperature;
	}

	friend std::ostream& operator << (std::ostream& o, const Edge& p)
	{	
		o << p.fromPhase << " --(T=" << p.temperature << ")--> " << p.toPhase;
		return o;
	}
};

struct Path
{
	std::vector<int> phases;
	// This will be in descending order, with one less element than phases. temperatures[i] will be the critical
	// temperature between phases[i] and phases[i+1].
	std::vector<double> temperatures;

	/*Path()
	{
		assert(false && "Should never be calling Path's default constructor!");
	}*/

	Path(int startPhase)
	{
		phases.push_back(startPhase);
		//phases.push_back(nextPhase);
		//temperatures.push_back(transitionTemperature);
	}

	// Copy constructor.
	Path(const Path& otherPath)
	{
		// These vectors will be copied to new vectors, which is good since we don't want changes to otherPath's
		// vectors to change this path's vectors.
		phases = otherPath.phases;
		temperatures = otherPath.temperatures;
	}

	void extend(Edge& edge)
	{
		assert(edge.fromPhase == phases.back() && "Attempted to extend the path using an edge that doesn't begin at the end of this path!");

		phases.push_back(edge.toPhase);
		temperatures.push_back(edge.temperature);
	}

	bool canUndergoTransition(Edge& edge)
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
			o << " --(T=" << p.temperatures[i-1] << ")--> " << p.phases[i];
		}

		return o;
	}
};

struct Vertex
{
	int phase;
	std::vector<Edge> edges;

	/*Vertex()
	{
		assert(false && "Should never be calling Vertex's default constructor!");
	}*/

	Vertex(int phase)
	{
		this->phase = phase;
	}

	void addEdge(Edge& edge)
	{
		assert(edge.fromPhase == phase && "Attempted to add an edge to a vertex that is not the start point of the edge!");
		edges.push_back(edge);
	}
};

//struct FrontierNode
//{
//	Edge& edge;
//	Path& path;
//
//	/*FrontierNode()
//	{
//		assert(false && "Should never be calling FrontierNode's default constructor!");
//	}*/
//
//	FrontierNode(const Edge& edge_, const Path& path_) :
//		edge(edge_), path(path_)
//	{
//		//
//	}
//	/*{
//		this->edge = edge_;
//		this->path = path_;
//	}*/
//};

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

	for(int i = 0; i < newPhase.X.size(); ++i)
	{
		for(int j = 0; j < reflectionIndices.size(); ++j)
		{
			if(!nonZero && abs(newPhase.X[i][reflectionIndices[j]]) > 1)
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
  *	Generates all combinations of the symmetryIndices list, stored in a list of lists. O(N*(2^N)), where N is the 
  * number of symmetric axes.
  */
std::vector<std::vector<int>> generateReflectionLists(const std::vector<int>& symmetryIndices)
{
	int N = symmetryIndices.size();
	int numSubLists = 1 << N;
	std::vector<std::vector<int>> reflectionLists;

	for(int i = 0; i < numSubLists; ++i)
	{
		std::vector<int> reflectionList = {};

		for(int j = 0; j < N; ++j)
		{
			if((i >> j) & 1)
			{
				reflectionList.push_back(symmetryIndices[j]);
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
			if(!nonZero && abs(newVacuum[reflectionIndices[i]]) > 1)
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
	const std::vector<int>& symmetryIndices,
	std::vector<PhaseTracer::Phase>& out_symmetrisedPhases,
	std::vector<PhaseTracer::Transition>& out_symmetrisedTransitions)
{
	// Generate all combinations of symmetryIndices.
	std::vector<std::vector<int>> reflectionLists = generateReflectionLists(symmetryIndices);

	// Grab a copy of the list of transitions.
	std::vector<PhaseTracer::Transition> tempTransitions = transitions;
	// We will use this filtered list henceforth.
	std::vector<PhaseTracer::Transition> filteredTransitions;

	// Remove any symmetric transitions, or rather keep only the transitions between the explictly stored phases.
	for(int i = 0; i < tempTransitions.size(); ++i)
	{
		if(tempTransitions[i].key == 0)
		{
			filteredTransitions.push_back(tempTransitions[i]);
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
	const std::vector<PhaseTracer::Transition>& transitions, Eigen::VectorXd& EWVEV)
{
	if(phases.size() == 0)
	{
		std::cout << "No phases identified." << std::endl;
		return PhaseStructureData(-1, -1, false);
	}

	if(transitions.size() == 0)
	{
		std::cout << "No transitions identified." << std::endl;
		return PhaseStructureData(-1, -1, false);
	}

	int EWVEVIndex = -1;
	int globalMinimumIndex = -1;

	double minPotentialValue = std::numeric_limits<double>::max();

	double maxT = std::numeric_limits<double>::min();
	double highTPhaseIndex = -1;


	// First identify the phase that corresponds to (vh, vs). Check all transitions that exist at T=0.
	for(int i = 0; i < phases.size(); ++i)
	{
		// At the start of this iteration, try to find if this phase is the high-T phase.
		// We check if it has the highest 'max T' (the highest temperature it is known to exist at), or if it is equal
		// to the current heighest max T, check if this is a lower energy phase than the other.
		if(phases[i].T.back() > maxT)
		{
			maxT = phases[i].T.back();
			highTPhaseIndex = i;
		}
		// Note: we are guaranteed to have highTPhaseIndex != -1 if we get here.
		else if(phases[i].T.back() == maxT && phases[i].V.back() < phases[highTPhaseIndex].V.back())
		{
			highTPhaseIndex = i;
		}
		
		// For the rest of this iteration we check the phases at T=0.	

		// Only check phases that exist at T=0. T.front() is the minimum temperature for the phase.
		if(phases[i].T.front() > 0)
		{
			continue;
		}

		Eigen::VectorXd zeroTPhase = phases[i].X.front();
		
		// Only worry about h > 0 phases, as the potential is symmetric in h and the EWVEV has h > 0.
		if(zeroTPhase[0] < -1.)
		{
			continue;
		}

		double zeroTPotential = phases[i].V.front();
		
		// Check if this is the deepest minimum (at T=0) found so far.
		if(zeroTPotential < minPotentialValue)
		{
			minPotentialValue = zeroTPotential;
			globalMinimumIndex = i;

			// If we have found the EWVEV already, then we know it's not the global minimum at T=0.
			if(EWVEVIndex != -1)
			{
				std::cout << "A deeper minimum than the EWVEV was found." << std::endl;
				return PhaseStructureData(EWVEVIndex, -1, false);
			}
		}

		if((zeroTPhase - EWVEV).squaredNorm() < 1.)
		{
			// If a deeper minimum has already been found and it's not another phase coinciding with the EWVEV, then
			// we know the EWVEV is not the global minimum at T=0.
			if(globalMinimumIndex < i && globalMinimumIndex != EWVEVIndex)
			{
				std::cout << "A deeper minimum than the EWVEV was found." << std::endl;
				return PhaseStructureData(EWVEVIndex, -1, false);
			}

			EWVEVIndex = i;
		}
	}

	if(EWVEVIndex == -1)
	{
		std::cout << "No EWVEV at T=0." << std::endl;
		return PhaseStructureData(-1, highTPhaseIndex, false);
	}

	assert(EWVEVIndex == globalMinimumIndex && "Logic to detect global minima not equal to the EWVEV is flawed!");
	
	return PhaseStructureData(EWVEVIndex, highTPhaseIndex, true);
}

std::vector<Path> getTransitionPaths(const std::vector<PhaseTracer::Phase>& phases,
	const std::vector<PhaseTracer::Transition>& transitions, const PhaseStructureData& phaseStructureData)
{	
	if(transitions.size() == 0)
	{
		std::vector<Path> paths;

		// If there are no transitions, but the high temperature phase becomes the EW VEV at T=0, then we have a valid
		// 'path' between the phases.
		if(phaseStructureData.highTPhaseIndex == phaseStructureData.EWVEVIndex)
		{
			paths.push_back(Path(phaseStructureData.highTPhaseIndex));
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
		Edge edge(transitions[i].false_phase.key, transitions[i].true_phase.key, transitions[i].TC);
		vertices[edge.fromPhase].addEdge(edge);
	}

	// Construct the initial frontier from all edges coming from the high temperature phase.
	int initialFrontierSize = vertices[phaseStructureData.highTPhaseIndex].edges.size();
	//std::vector<FrontierNode> frontier(initialFrontierSize);
	std::vector<FrontierNode> frontier;
	int maxTemperature = phases[phaseStructureData.highTPhaseIndex].T.back();
	//std::vector<Path> paths(initialFrontierSize);
	std::vector<Path> paths;

	for(int i = 0; i < initialFrontierSize; ++i)
	{
		//Edge& edge = vertices[phaseStructureData.highTPhaseIndex].edges[i];
		//Path path(edge.fromPhase);
		//paths[i] = path;
		Path path(phaseStructureData.highTPhaseIndex);
		paths.push_back(path);
		//frontier[i] = FrontierNode(edge, path);
		//FrontierNode frontierNode(edge, path);
		FrontierNode frontierNode(phaseStructureData.highTPhaseIndex, i, i);
		frontier.push_back(frontierNode);
	}
	
	// Use a breadth-first search to find all edges from the high temperature phase to the EW VEV, using a frontier
	// approach. When some the current end point of some path has multiple valid edges from it, we need to create new
	// paths that branch from this point, but are identical up to this current point.
	while(frontier.size() > 0)
	{
		std::vector<FrontierNode> newFrontier;

		for(int i = 0; i < frontier.size(); ++i)
		{
			Edge& edge = vertices[frontier[i].vertexIndex].edges[frontier[i].edgeIndex];
			int vertexIndex = edge.toPhase;
			Path& path = paths[frontier[i].pathIndex];

			path.extend(edge);

			// This has been removed since we do care about extending the path further if there is a valid transition
			// away from the EW VEV. Of course, we need to relax conditions in extractPhaseStructureData for this to be
			// considered anyway.
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

	// Check that all paths end at the EW VEV.
	/*for(int i = 0; i < paths.size(); ++i)
	{
		bool valid = paths[i].phases.back() == phaseStructureData.EWVEVIndex;
		std::cout << (valid ? "[Valid]  " : "[Invalid]") << " Path " << i+1  << ": " << paths[i] << std::endl;
		//assert(paths[i].phases.back() == phaseStructureData.EWVEVIndex && "A path did not end at the EW VEV!");
	}*/

	return paths;
}

} // namespace TransitionGraph

#endif
