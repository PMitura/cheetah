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
        virtual Points2D solve(const PointsND& inputSet);

    protected:
        /**
         * Internal represenation of algorithm for finding the convex hull
         */
        virtual Points2D algorithm(const data_t& input) = 0;
};

}
