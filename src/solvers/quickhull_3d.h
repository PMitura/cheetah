#pragma once

#include <vector>

#include "lib/geometry.h"
#include "lib/structures.h"
#include "solvers/solver_3d.h"

namespace ch
{

/**
 * Quickhull modified for convex hulls in 3D
 */
class Quickhull3D : public Solver3D
{
    public:
        Quickhull3D();
        /**
         * Find convex hull of input data set.
         *
         * @param inputSet Input data set
         * @return Convex hull of input data set.
         */
        Polyhedron& solve(const Points3D& input, Polyhedron& output);

    private:
        Polyhedron& solveSequential(const Points3D& input, Polyhedron& output);

        /** finds initial tetrahedron */
        Polyhedron findInitial(const Points3D& input);


};

}
