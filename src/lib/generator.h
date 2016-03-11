#pragma once

#include <algorithm>
#include <cstdlib>
#include <ctime>

#include "lib/structures.h"
#include "lib/geometry.h"

namespace ch
{

class Generator2D
{
    public:
        /**
         * Generates random data with specified parameters
         * 
         * Hull points are spread uniformly on circle specified by radius,
         * interior points are generated as convex combination of hull points
         * or their subset
         *
         * @param n      Total number of points in dataset
         * @param h      Number of points on hull
         * @param radius Range of coordinates points can have
         * @param points Output set of points
         * @return True on success, false on invalid parameters (e.g. n < h)
         */
        bool genUniformCircle(long long n, long long h,
                              double radius, Points2D& points);

        /**
         * Same as genUniformCircle, only the hull points are spread randomly
         */
        bool genRandomCircle(long long n, long long h,
                              double radius, Points2D& points);

    private:
        /**
         * Generates points in linear combination of given points
         *
         * @param n      Number of points to generate
         * @param of     Number of points in combination
         * @param points Contains points in combination, serves as output
         */
        void generateInCombination(long long n, int of, Points2D& points);

        /**
         * Triangulates given convex polygon, then generates points inside
         * these triangles unifromly.
         *
         * Slower than combination, but truly uniform
         *
         * @param n      Number of points to generate
         * @param of     Number of points in polygon
         * @param points Contains points in polygon, serves as output
         */
        void generateTriangulated(long long n, int of, Points2D& points);
};

}
