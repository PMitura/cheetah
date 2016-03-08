#pragma once

#include <iostream>
#include <vector>

namespace ch
{

typedef std::vector<std::vector<double>> data_t;

/**
 * General class for input or output point set
 */
class PointsND
{
    public:
        /**
         * Constructor.
         *
         * @param dim Dimension of point set
         */
        PointsND(unsigned int dim);

        /**
         * Add a point to the set
         */
        bool add(std::vector<double> point);

        /** Get number of points in set */
        inline int getSize() const { return data_.size(); }

        /** Data getter */
        inline const data_t& getData() const { return data_; }

        /** Dimension getter */
        inline unsigned int getDimension() const { return dimension_; }

    protected:
        /** Internal representation of points */
        data_t data_;

        /** Dimension of point set */
        unsigned int dimension_;
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
