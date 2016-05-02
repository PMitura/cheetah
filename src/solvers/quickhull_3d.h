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
        /** Beginning point for sequential solving */
        Polyhedron& solveSequential(const Points3D& input, Polyhedron& output);

        /** 
         * Finds initial tetrahedron.
         *
         * @param input all input pts
         * @return tetrahedron, or empty polyhedron in case of degenerate input
         */
        Polyhedron findInitial(const data_t& input);

        struct QVertex;
        struct QHalfEdge;
        struct QFace;

        struct QVertex {
            QVertex(point_t pt)
                : next_(NULL), prev_(NULL), edge_(NULL), crds_(pt) {}
            QVertex() : QVertex({0, 0, 0}) {}
            QVertex * next_, * prev_;
            QHalfEdge * edge_;
            point_t crds_;
            double operator[](unsigned x) {
                return crds_[x];
            }
        };

        struct QHalfEdge {
            QHalfEdge() : tail_(NULL), twin_(NULL), prev_(NULL), next_(NULL),
                  face_(NULL) {}
            QVertex * tail_;
            QHalfEdge * twin_, * prev_, * next_;
            QFace * face_;
        };

        struct QFace {
            QFace() : next_(NULL), prev_(NULL), edge_(NULL) {}
            QFace * next_, * prev_;
            QHalfEdge * edge_;
            std::vector<unsigned> conflicts_;
        };
};

}
