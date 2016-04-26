#pragma once

#include <vector>
#include <cstdlib>
#include <ctime>
#include <cfloat>

#include "lib/geometry.h"
#include "lib/structures.h"
#include "solvers/solver_3d.h"

namespace ch
{

/**
 * Gift wrapping style solver for 3d convex hulls
 */
class JarvisScan3D : public Solver3D
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
        /** Returns 1 or -1 at random */
        double randomOne();
        std::pair<unsigned, unsigned> findInitial(const data_t& input);
};

}
