#pragma once

#include <iomanip>
#include <vector>
#include <stack>
#include <list>
#include <omp.h>

#include "lib/structures.h"
#include "lib/geometry.h"
#include "solvers/solver_2d.h"

namespace ch
{

class Quickhull2D : public Solver2D
{
    public:
        Quickhull2D();
        Points2D& solve(const Points2D& input, Points2D& output);

        enum Variant {NAIVE, FORWARD, PRECOMP, PARA};
        Quickhull2D(Variant v);
        Quickhull2D(Variant v, int threshold);

    private:
        void recNaive(point_t& a, point_t& b, data_t& plane);
        Points2D& solveNaive(const Points2D& input, Points2D& output);

        void recSequential(point_t& a, point_t& b, point_t& c, data_t& plane);
        void recSplit(point_t& a, point_t& b, point_t& c, data_t& plane,
                      bool upper);

        Points2D& solveSequential(const Points2D& input, Points2D& output);

        void recPrecomp(point_t& a, point_t& b, point_t& c, data_t& plane);
        Points2D& solvePrecomp(const Points2D& input, Points2D& output);

        void recForwarded(const point_t& a, const point_t& b, const point_t& c, 
                         std::vector<unsigned> plane, unsigned planeSize);
        Points2D& solveForwarded(const Points2D& input, Points2D& output);

        void recParallel(const point_t& a, const point_t& b, const point_t& c, 
                         std::vector<unsigned> plane, unsigned planeSize,
                         std::list<point_t>& onHull);
        Points2D& solveParallel(const Points2D& input, Points2D& output);

        Points2D& solveIterative(const Points2D& input, Points2D& output);

        std::pair<point_t, point_t> minMaxX(const data_t& points);
        std::pair<point_t, point_t> farthestPoints(const data_t& points);

        point_t planeFarthestCross(point_t& a, point_t& b,
                                   const data_t& plane);

        point_t planeFarthestDist(point_t& a, point_t& b,
                                   const data_t& plane);

        void divideToPlanes(const data_t& input,
                            point_t& pivotLeft, point_t& pivotRight,
                            data_t& topPlane, data_t& botPlane);
        void divideToPlanesPara(const data_t& input,
                                point_t& pivotLeft, point_t& pivotRight,
                                data_t& topPlane, data_t& botPlane);

        unsigned int parallelThreshold_;
        const data_t* globIn_;
        Points2D* globOut_;

        /** local epsilon value */
        double EPS_LOC;

        Variant variant_;
};

}
