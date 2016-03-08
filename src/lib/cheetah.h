#pragma once

#include <iostream>
#include <vector>

namespace ch
{

/**
 * General class for input or output point set
 */
class PointsND
{
    public:
        /**
         * Constructor 
         */
        PointsND(int dim);

        /** Dimension getter */
        inline int getDimension() { return dimension_; }

    protected:
        /** Internal representation of points */
        std::vector<std::vector<double>> data_;

        /** Dimension of point set */
        int dimension_;
};

/**
 *
 */
class Points2D : public PointsND
{
    public:
        Points2D();
};

}
