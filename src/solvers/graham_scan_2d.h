#pragma once

#include "solvers/solver_2d.h"

namespace ch
{

class GrahamScan2D : public Solver2D
{
    public:
        points2D solve(const points2D& inputSet);
};

}

