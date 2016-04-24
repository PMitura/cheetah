#pragma once

#include "lib/structures.h"

namespace ch
{

/**
 * Abstract class serving as a template for 3d convex hull solver.
 */
class Solver3D
{
    public:
        virtual ~Solver3D() {}

        /**
         * Find convex hull of input data set.
         *
         * @param inputSet Input data set
         * @return Convex hull of input data set.
         */
        virtual Polyhedron& solve(const Points3D& input, Polyhedron& output) = 0;

        inline std::string getName() { return name_; }

    protected:
        /** name of solver */
        std::string name_;
};

}
