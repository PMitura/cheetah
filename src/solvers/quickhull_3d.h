#pragma once

#include <vector>

#include "lib/geometry.h"
#include "lib/half_edge.h"
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
        /** Beginning point for sequential solving */
        Polyhedron& solveSequential(const Points3D& input, Polyhedron& output);

        /** 
         * Finds initial tetrahedron.
         *
         * @param input all input pts
         * @return tetrahedron, or empty polyhedron in case of degenerate input
         */
        void findInitial(const data_t& input);

        /**
         * Global list of found faces in half-edge mesh structure
         */
        std::vector<QFace> faces_;
};

}
