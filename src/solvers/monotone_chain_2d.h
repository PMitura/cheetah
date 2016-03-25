#pragma once

#include <vector>
#include <algorithm>
#include <omp.h>
#include <parallel/algorithm>
#include <parallel/settings.h>

#include "solvers/solver_2d.h"
#include "lib/structures.h"
#include "lib/geometry.h"

namespace ch
{

class MonotoneChain2D : public Solver2D
{
    public:
        MonotoneChain2D();

        Points2D& solve(const Points2D& input, Points2D& output);

    private:
        Points2D& solveSequential(const Points2D& input, Points2D& output);
        Points2D& solveParallel(const Points2D& input, Points2D& output);

        unsigned scanLower(const data_t& input, unsigned* lower);
        unsigned scanUpper(const data_t& input, unsigned* upper);

        void sortPtsDirect(const data_t& input);
        void sortPtsCache(const data_t& input);
        void sortPtsParallel(const data_t& input);

        std::vector<unsigned> order_;
        std::vector<double> xar_, yar_;

        struct PointCmpDirect {
            PointCmpDirect(const data_t& d)
                : data_(d) {}
            bool operator()(const unsigned& a, const unsigned& b);
            const data_t& data_;
        };

        struct PointCmpCache {
            PointCmpCache(const MonotoneChain2D& p, const data_t& d)
                : part_(p), data_(d) {}
            bool operator()(const unsigned& a, const unsigned& b);
            const MonotoneChain2D& part_;
            const data_t& data_;
        };
};

}
