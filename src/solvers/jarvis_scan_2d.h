#pragma once

#include <omp.h>

#include "lib/structures.h"
#include "lib/geometry.h"
#include "solvers/solver_2d.h"

namespace ch
{

class JarvisScan2D : public Solver2D
{
    public:
        JarvisScan2D();

        Points2D& solve(const Points2D& input, Points2D& output);

        /** Algorithm version switcher */
        enum Variant {CROSS, POLAR, PARA, PARA_INT, PARA_DOUBLE};
        JarvisScan2D(Variant v);
        inline void setVariant(Variant v) { variant_ = v; }

    private:
        Points2D& solveCross(const Points2D& input, Points2D& output);
        Points2D& solvePolar(const Points2D& input, Points2D& output);
        Points2D& solvePara(const Points2D& input, Points2D& output);

        void scan(const data_t& input, data_t& output,
                  unsigned beginIdx, unsigned endIdx);
        void scanPara(const data_t& input, data_t& output,
                  unsigned beginIdx, unsigned endIdx);

        Variant variant_;
};

}
