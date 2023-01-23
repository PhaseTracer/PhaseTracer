#ifndef TRANSITIONGRAPH_TRANSITION_GRAPH_UTIL_HPP_INCLUDED
#define TRANSITIONGRAPH_TRANSITION_GRAPH_UTIL_HPP_INCLUDED

//#include "logger.hpp"
#include "transition_finder.hpp"
#include "phase_finder.hpp"

#include <assert.h>
#include <stdlib.h>
#include <limits>

// Forward declarations.
namespace PhaseTracer
{
	struct Transition;
	class TransitionFinder;
}

namespace TransitionGraph
{

/**
  * Stores the boundary conditions of the phase history (i.e. the phases that exist at the high temperature, and the
  * lowest energy phases at the low temperature (T=0)), and whether the model is valid at the low temperature (i.e.
  * whether the lowest energy phase is in the expected field configuration).
  *
  * This struct is constructed in -- and returned from -- TransitionGraph::extractPhaseStructureData.
  */
struct PhaseStructureData
{
	std::vector<int> lowTPhaseIndices;
	std::vector<int> highTPhaseIndices;
	std::vector<bool> isLowTemperaturePhase;
	bool validAtZeroT;

	PhaseStructureData(std::vector<int> lowTPhaseIndices, std::vector<int> highTPhaseIndices, bool validAtZeroT,
		std::vector<bool> isLowTemperaturePhase)
	{
		this->lowTPhaseIndices = lowTPhaseIndices;
		this->highTPhaseIndices = highTPhaseIndices;
		this->validAtZeroT = validAtZeroT;
		this->isLowTemperaturePhase = isLowTemperaturePhase;
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

struct TransitionEdge
{
	int transitionIndex;
	bool subcritical;
	double temperature;

	TransitionEdge(int transitionIndex, bool subcritical, double temperature)
	{
		this->transitionIndex = transitionIndex;
		this->subcritical = subcritical;
		this->temperature = temperature;
	}

	friend std::ostream& operator << (std::ostream& o, const TransitionEdge& te)
	{
		o << "{i: " << te.transitionIndex << ", sc: " << te.subcritical << ", T: " << te.temperature << "}";

		return o;
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
		//LOG(debug) << "Adding edge " << edge << " to vertex " << phase;
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

struct Path
{
	std::vector<int> phases;
	// This will be in descending order, with one less element than phases. temperatures[i] will be the critical
	// temperature between phases[i] and phases[i+1].
	//std::vector<double> temperatures;
	// Similar to temperatures, this will be the transition indices in descending order of the corresponding
	// temperatures.
	std::vector<TransitionEdge> transitions;

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

	void extend(const Edge& edge, bool subcritical, double temperature)
	{
		assert(edge.fromPhase == phases.back() && "Attempted to extend the path using an edge that doesn't begin at the end of this path!");

		phases.push_back(edge.toPhase);
		//temperatures.push_back(edge.temperature);
		//transitions.push_back(edge.transition);
		transitions.push_back({edge.transition, subcritical, temperature});
	}

	//bool canUndergoTransition(const Edge& edge)
	bool canUndergoTransition(const std::vector<Vertex>& vertices, int vertexIndex, int edgeIndex)
	{
		// TODO: need to make sure this function isn't used for checking higher temperature subcritical transitions from
		// the new toPhase. If it is, this condition needs to be changed.
		//return edge.temperature <= temperatures.back();
		//return edge.temperature <= transitions.back().temperature;

		// The phase can undergo this transition if either the transition temperature is below this path's current
		// temperature, or if there does not exist the reverse transition below it's temperature. Additionally, it should
		// not return the system to a phase already existing in this path but at a lower temperature, as this new
		// transition at a higher temperature .

		assert(vertexIndex == phases.back() && "Checked whether a transition was possible from a phase that isn't the end of this path!");

		const Vertex& vertex = vertices[vertexIndex];
		const Edge& transition = vertex.edges[edgeIndex];

		int fromPhase = vertex.edges[edgeIndex].fromPhase;
		int toPhase = vertex.edges[edgeIndex].toPhase;
		double currentTemperature = transitions.back().temperature;

		// If this transition occurs at a lower temperature than the path's current temperature, then the transition
		// can occur.
		if(transition.temperature <= currentTemperature)
		{
			return true;
		}

		// If this transition occurs at a higher temperature than the path's current temperature, then there are
		// conditions that must be satisfied for the transition to be possible.


		// Check if there is a reverse transition between the current temperature and the transition temperature. If
		// there is, this transition cannot occur.
		const Vertex& toVertex = vertices[transition.toPhase];

		for(int i = 0; i < toVertex.edges.size(); ++i)
		{
			if(toVertex.edges[i].temperature < currentTemperature)
			{
				continue;
			}

			if(toVertex.edges[i].temperature > transition.temperature)
			{
				break;
			}

			if(toVertex.edges[i].toPhase == vertex.phase)
			{
				return false;
			}
		}

		// Now we need to check if this transition returns the path to a phase already existing in the path but at a 
		// lower temperature. If we return to a phase already existing in the path, but return to it at a higher
		// temperature, then we should discard this transition. We have already handled that phase at the higher
		// temperature, and doing so again would result in cycles.
		
		// Since the path is in order of monotonically decreasing temperature, if we find the phase in the path at all
		// then we know it must be at a higher temperature.
		for(int i = 0; i < phases.size(); ++i)
		{
			if(phases[i] == toPhase)
			{
				return false;
			}
		}
		
		// Otherwise, this transition can be added.
		return true;
	}

	// Returns the current temperature of the path, given by the temperature at which the last transition occurred.
	// This is used for determining the temperature at which a newly added transition will occur. For instance, if the
	// newly added transition has a critical temperature above the path's current temperature, then it will be treated as
	// a subcritical transition with a transition temperature equal to the path's current temperature.
	double getCurrentTemperature()
	{
		return transitions.size() > 0 ? transitions.back().temperature : std::numeric_limits<double>::max();
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
			//o << " --(" << p.transitions[i-1].tr << ", T=" << p.temperatures[i-1] << ")--> " << p.phases[i];
			o << " --" << p.transitions[i-1] << "--> " << p.phases[i];
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

		//transitionIndices += std::to_string(transitions[0]);
		transitionIndices += std::to_string(transitions[0].transitionIndex);

		for(int i = 1; i < transitions.size(); ++i)
		{
			//transitionIndices += " " + std::to_string(transitions[i]);
			transitionIndices += " " + std::to_string(transitions[i].transitionIndex);
		}

		return transitionIndices;
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
	reflectionIndices, int key);

PhaseTracer::Transition constructSymmetricPartnerTransition(const PhaseTracer::Transition& transition, const 
	PhaseTracer::Phase& false_phase, const PhaseTracer::Phase& true_phase, const Eigen::VectorXd& false_vacuum,
	const Eigen::VectorXd& true_vacuum);

/**
  *	Generates all combinations of the symmetryIndices list, stored in a list of lists. O(n*N*(2^N)), where N is the 
  * number of discrete symmetries and n is the number of fields.
  */
std::vector<std::vector<int>> generateReflectionLists(const std::vector<std::vector<int>>& symmetryIndices);

/** Returns the symmetric vacuum using the axes of reflection defined through reflectionIndices. */
Eigen::VectorXd getReflectedVacuum(const Eigen::VectorXd& vacuum, const std::vector<int>& reflectionIndices,
	bool checkRedundancy);

/** The comparison function for sorting phases by their keys. */
bool comparePhases(const PhaseTracer::Phase& a, const PhaseTracer::Phase& b);

/** The comparison function for sorting transitions by their false vacuum keys. */
bool compareTransitions(const PhaseTracer::Transition& a, const PhaseTracer::Transition& b);

// TODO: fix excessive copying?
void extractExplicitSymmetricPhasesAndTransitions(
	const std::vector<PhaseTracer::Phase>& phases,
	const std::vector<PhaseTracer::Transition>& transitions,
	const std::vector<std::vector<int>>& symmetryIndices,
	std::vector<PhaseTracer::Phase>& out_symmetrisedPhases,
	std::vector<PhaseTracer::Transition>& out_symmetrisedTransitions);

PhaseStructureData extractPhaseStructureData(const std::vector<PhaseTracer::Phase>& phases,
	const std::vector<PhaseTracer::Transition>& transitions, const std::vector<Eigen::VectorXd>& expectedLowTPhases,
	double Tmax, bool knownHighTPhase);

std::vector<Path> getTransitionPathsFromHighTPhase(const std::vector<Vertex>& vertices, int highTPhaseIndex);

std::vector<Path> getTransitionPaths(const std::vector<PhaseTracer::Phase>& phases,
	const std::vector<PhaseTracer::Transition>& transitions, const PhaseStructureData& phaseStructureData);

/**
  * knownHighTPhase=true means we know the global minimum at T=Tmax is the phase the Universe is in at T=Tmax.
  * knownHighTPhase=false means we don't know which of the phases at T=Tmax the Universe is in at T=Tmax.
  * This boolean variable is equivalent to whether we have sampled the potential at high enough temperatures.
  */
std::vector<Path> getPhaseHistory(const PhaseTracer::PhaseFinder& pf, const
	PhaseTracer::TransitionFinder& tf, const EffectivePotential::Potential& model, bool knownHighTPhase);

} // namespace TransitionGraph

#endif
