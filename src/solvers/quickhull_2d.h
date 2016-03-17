#pragma once

#include "lib/structures.h"
#include "lib/geometry.h"
#include "solvers/solver_2d.h"

namespace ch {

class Quickhull2D : public Solver2D
{
    public:
        Quickhull2D();
        Points2D& solve(const Points2D& input, Points2D& output);

    private:
        void recNaive(point_t& a, point_t& b, data_t& plane);
        Points2D& solveNaive(const Points2D& input, Points2D& output);

        Points2D* globOut_;
};

}