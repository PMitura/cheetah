#pragma once

#include "lib/structures.h"
#include "lib/geometry.h"
#include "solvers/solver_2d.h"

namespace ch {

class Quickhull2D : public Solver2D
{
    public:
        Points2D& solve(const Points2D& input, Points2D& output);

    private:
        Points2D& solveNaive(const Points2D& input, Points2D& output);
};

}
