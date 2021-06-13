#pragma once

#include "Matrix2D.h"
#include "Reachability.h"

#include <vector>


class MarkovChain {
public:
	MarkovChain() {}

	Mesh p_t_;
	Mesh p_;
};

MarkovChain CalculateMarkovChainProbabilities(const std::vector<Reachability>& reachability, const Scenario& scenario, const Mesh& previous_state);

