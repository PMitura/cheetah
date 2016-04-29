#pragma once

#include <vector>

#include "lib/geometry.h"
#include "lib/structures.h"
#include "solvers/solver_3d.h"

namespace ch
{

/**
 * Quickhull modified for convex hulls in 3d
 */
class Quickhull3D : public Solver3D
{
    public:
        /**
         * Find convex hull of input data set.
         *
         * @param inputSet Input data set
         * @return Convex hull of input data set.
         */
        Polyhedron& solve(const Points3D& input, Polyhedron& output);

    private:
        Polyhedron& solveNaive(const Points3D& input, Polyhedron& output);
};

}
