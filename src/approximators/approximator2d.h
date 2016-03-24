#pragma once

#include "lib/structures.h"

namespace ch
{

/**
 * Abstract class serving as a template for 2d convex hull solver.
 */
class Appproximator2D
{
    public:
        virtual ~Appproximator2D() {}

        /**
         * Find convex hull of input data set.
         *
         * @param inputSet Input data set
         * @return Convex hull of input data set.
         */
        virtual Points2D& approximate(const Points2D& input,
                                      Points2D& output) = 0;

        inline std::string getName() { return name_; }

    protected:
        /** name of solver */
        std::string name_;
};

}
