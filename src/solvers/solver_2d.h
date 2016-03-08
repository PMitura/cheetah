#pragma once

#include <vector>
#include <lib/cheetah.h>

namespace ch
{

/**
 * Abstract class serving as a template for 2d convex hull solver.
 */
class Solver2D
{
    public:
        /**
         * Find convex hull of input data set.
         *
         * @param inputSet Input data set
         * @return Convex hull of input data set.
         */
        virtual Points2D& solve(const Points2D& input, Points2D& output) = 0;
};

}
