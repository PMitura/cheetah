#pragma once

#include "lib/structures.h"
#include "lib/geometry.h"
#include "solvers/solver_2d.h"

namespace ch {

class JarvisScan2D : public Solver2D
{
    protected:
        Points2D& solve(const Points2D& input, Points2D& output);
};

}
