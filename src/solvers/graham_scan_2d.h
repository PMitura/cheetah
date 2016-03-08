#pragma once

#include "solvers/solver_2d.h"

namespace ch
{

class GrahamScan2D : public Solver2D
{
    public:
        Points2D solve(const Points2D& inputSet);
};

}

