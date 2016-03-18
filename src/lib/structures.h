#pragma once

#include <iostream>
#include <vector>
#include <utility>

namespace ch
{

typedef std::vector<double> point_t;
typedef std::vector<point_t> data_t;
typedef std::pair<double, double> point2d_t;

// debug macro, should not appear in relase
#define D(X) std::cout<<"  "<<#X": "<<X<<std::endl;
#define R(X) std::cout<<"  "<<X<<std::endl;

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
        bool add(point_t point);

        /** Remove all points from set */
        inline void clear() { data_.clear(); }

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
