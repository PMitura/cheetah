#pragma once

#include <vector>
#include <stack>

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

    private:
        void recNaive(point_t& a, point_t& b, data_t& plane);
        Points2D& solveNaive(const Points2D& input, Points2D& output);
        Points2D& solveIterative(const Points2D& input, Points2D& output);

        std::pair<point_t, point_t> minMaxX(const data_t& points);
        std::pair<point_t, point_t> farthestPoints(const data_t& points);

        point_t planeFarthestPoint(point_t& a, point_t& b,
                                   const data_t& plane);

        void divideToPlanes(const data_t& input,
                            point_t& pivotLeft, point_t& pivotRight,
                            data_t& topPlane, data_t& botPlane);

        Points2D* globOut_;
};

}
