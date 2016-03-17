#pragma once

#include <iostream>
#include <vector>
#include <ctime>
#include <omp.h>

#include "lib/generator.h"
#include "lib/structures.h"
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

        bool runGeneratedTest(int n, int h, double span, Solver2D& solver);

    private:

};

}
