#pragma once

#include "lib/structures.h"

namespace ch
{

/**
 * Abstract class serving as a template for 2d convex hull approximator.
 */
class Approximator2D
{
    public:
        virtual ~Approximator2D() {}

        /**
         * Find convex hull approximation of input data set.
         *
         * @param inputSet Input data set
         * @return convex hull approximation of input data set.
         */
        virtual Points2D& approximate(const Points2D& input,
                                      Points2D& output) = 0;

        inline std::string getName() { return name_; }

        /** Returns maximum reachable number of points on output hull */
        virtual int maxReachable() const = 0;

    protected:
        /** name of approximator */
        std::string name_;
};

}
