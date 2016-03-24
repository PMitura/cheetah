#pragma once

#include <vector>
#include <algorithm>
#include <omp.h>
#include <parallel/algorithm>
#include <parallel/settings.h>

#include "solvers/solver_2d.h"
#include "lib/structures.h"
#include "lib/geometry.h"

namespace ch
{

class MonotoneChain2D : public Solver2D
{
    public:
        MonotoneChain2D();

        Points2D& solve(const Points2D& input, Points2D& output);

    private:
        Points2D& solveSequential(const Points2D& input, Points2D& output);
        Points2D& solveParallel(const Points2D& input, Points2D& output);

        unsigned scanLower(const data_t& input, unsigned* lower);
        unsigned scanUpper(const data_t& input, unsigned* upper);

        struct PointCmp {
            bool operator()(const point_t& a, const point_t& b);
        };
};

}
