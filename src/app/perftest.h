#pragma once

#include "solvers/solver_2d.h"
#include "solvers/jarvis_scan_2d.h"
#include "solvers/graham_scan_2d.h"
#include "solvers/quickhull_2d.h"

namespace ch
{

/**
 * Wrapper for running performance tests of implemented algorithms
 */
class PerfTest
{
    public:
        /** runs all scheduled tests */
        void runAllTests();

        bool runGeneratedTest(int n, int h, int dim);

    private:

};

}
