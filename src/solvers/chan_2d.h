#pragma once

#include <vector>
#include <cfloat>
#include <omp.h>

#include "lib/structures.h"
#include "solvers/graham_scan_2d.h"
#include "solvers/solver_2d.h"

namespace ch
{

class Chan2D : public Solver2D
{
    public:
        Chan2D();
        Points2D& solve(const Points2D& input, Points2D& output);

    private:
        Points2D& solveNaive(const Points2D& input, Points2D& output);
        void findHulls(const Points2D& input, std::vector<Points2D>& hulls,
                       unsigned hullSize);

        /**
         * Finds point with minimal x coordinate, along with hull it lies on
         *
         * @param hulls vector of hulls
         * @param minPt point with min x
         * @return index of hull min point lies on
         */
        unsigned findMinHull(std::vector<Points2D&> hulls, point_t& minPt);
};

}
