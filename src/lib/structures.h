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
#define D(X) //std::cout<<"  "<<#X": "<<X<<std::endl;
#define R(X) //std::cout<<"  "<<X<<std::endl;

const unsigned MAX_NUM_THREADS = 24;

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
        inline unsigned getSize() const { return data_.size(); }

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

class Points2D : public PointsND
{
    public:
        Points2D();
};

class Points3D : public PointsND
{
    public:
        Points3D();
};

class Polyhedron
{
    public:
        inline void addFace(Points3D& face) { faces_.push_back(face); }
        inline void popFace() { faces_.pop_back(); }
        inline const std::vector<Points3D>& getFaces() { return faces_; }
        inline unsigned getSize() { return faces_.size(); }

    private:
        std::vector<Points3D> faces_;
};


}
