#pragma once

#include <vector>

#include "solvers/solver_2d.h"
#include "lib/geometry.h"

namespace ch
{

class GrahamScan2D : public Solver2D
{
    public:
        GrahamScan2D();
        ~GrahamScan2D();

        Points2D& solve(const Points2D& input, Points2D& output);

    private:
        /** finds point with minimum Y in given set */
        int findMinY(const data_t& points);

        /** Precomputes polar angles of all points, with respect to minY */
        void computeAngles(int pivotIdx);

        /** point indexes sorted by polar angle */
        std::vector<int> * order_;
        /** precomputed polar angles */
        std::vector<double> * polar_;
};

}

