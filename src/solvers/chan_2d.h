#pragma once

#include <vector>
#include <cfloat>
#include <omp.h>

#include "lib/structures.h"
#include "solvers/jarvis_scan_2d.h"
#include "solvers/graham_scan_2d.h"
#include "solvers/monotone_chain_2d.h"
#include "solvers/quickhull_2d.h"
#include "solvers/solver_2d.h"

namespace ch
{

class Chan2D : public Solver2D
{
    public:
        Chan2D();
        Points2D& solve(const Points2D& input, Points2D& output);

        enum Variant {JARVIS, GRAHAM, QUICK, COMBO};
        Chan2D(Variant v);

    private:
        inline unsigned ppow(unsigned x) { return 1U << (1U << x); }
        Points2D& solveNaive(const Points2D& input, Points2D& output);
        void findHulls(const Points2D& input, std::vector<Points2D>& hulls,
                       unsigned step);

        /**
         * Finds point with minimal x coordinate, along with hull it lies on
         *
         * @param hulls vector of hulls
         * @param minPt point with min x
         * @return index of hull min point lies on
         */
        unsigned findMinHull(std::vector<Points2D>& hulls, unsigned& minPt);

        /**
         * Finds next point on overall hull
         *
         * @param hulls list of subhulls
         * @param pt pair of point index and hull index of current point
         * @return pair of point index and hull index
         */
        std::pair<unsigned, unsigned> findNext(std::vector<Points2D>& hulls,
                std::pair<unsigned, unsigned> curr);

        /**
         * Finds point on hull touched by left tangent line from point p
         *
         * Complexity O(n log n) using binary search
         *
         * @param hull ordered convex subhull we are searching on
         * @param p specified point of tangent line
         * @return id of found point on its subhull
         */
        unsigned findTangent(const Points2D& hull, point_t& p);

        Variant variant_;
        Solver2D* solver_;
        bool comboFlag_;
};

}
