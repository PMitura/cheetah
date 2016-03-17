#pragma once

#include <vector>
#include <algorithm>

#include "solvers/solver_2d.h"
#include "lib/geometry.h"

namespace ch
{

class GrahamScan2D : public Solver2D
{
    public:
        GrahamScan2D();

        Points2D& solve(const Points2D& input, Points2D& output);

    private:
        /** finds point with minimum Y in given set */
        int findMinY(const data_t& points);

        /** Precomputes polar angles of all points, with respect to minY */
        void computeAngles(const data_t& points);

        /** point indexes sorted by polar angle */
        std::vector<unsigned> order_;
        /** precomputed polar angles */
        std::vector<double> polar_;
        /** pivot point index */
        int pivot_;

        struct AngleCmp {
            AngleCmp(const GrahamScan2D& p, const data_t& d)
                : part_(p), data_(d) {}
            bool operator()(const int& a, const int& b);
            const GrahamScan2D& part_;
            const data_t& data_;
        };
};

}

