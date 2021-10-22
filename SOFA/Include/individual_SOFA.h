#pragma once
#include "Individual.h"

namespace GlobalOptimization {
class IndividualSOFA : public Individual {
public:
 double probability;
 double cumulative_probability;
 IndividualSOFA& operator=(IndividualSOFA& ind);
};
}

