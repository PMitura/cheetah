#pragma once

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <ctime>
#include <omp.h>

#include "cheetah/core.h"

#include "approximators/approximator2d.h"
#include "approximators/bfp2d.h"
#include "lib/generator.h"
#include "lib/structures.h"
#include "solvers/solver_2d.h"
#include "solvers/solver_3d.h"
#include "solvers/jarvis_scan_2d.h"
#include "solvers/jarvis_scan_3d.h"
#include "solvers/graham_scan_2d.h"
#include "solvers/monotone_chain_2d.h"
#include "solvers/quickhull_2d.h"
#include "solvers/chan_2d.h"

namespace ch
{

/**
 * Wrapper for running performance tests of implemented algorithms
 */
class PerfTest
{
    public:
        /** Runs all scheduled tests */
        void runAllTests();

        /** Runs one test with specified parameters */
        double runTest(int n, int h, int span, int runs,
                int cores, SolverType type);

    private:
        /** One instance of test run */
        struct Instance {
            int n, h, runs;
            double span;
            int cores;
        };

        /**
         * Runs series of small instance tests with various parameters.
         *
         * Tests are run multiple times with same parameters to reduce
         * inacurracies in measurement.
         *
         * @param solvers list of tested solvers
         */
        void smallTests(std::vector<Solver2D*> solvers);

        /**
         * Runs series of tests on big inputs with various parameters.
         *
         * @param solvers list of tested solvers
         */
        void bigTests(std::vector<Solver2D*> solvers);

        /**
         * Runs series of tests to find edges for various input sizes.
         */
        void testEdges();

        /**
         * Finds point where QH beats JS, and GS beats QH. If that happens.
         *
         * @param setSize number of points in input set
         */
        void findEdge(int setSize);

        /**
         * Runs single instance of test, on group of solvers.
         *
         * @param inst instance parameters
         * @param solvers solvers to use on specified instance
         */
        void runTestInstance(Instance& inst, std::vector<Solver2D*> solvers);

        /**
         * Runs, and measures single problem instance.
         *
         * Time to generate data is not counted.
         *
         * @param h target number of points on hull
         * @param input set of points to solve
         * @param solver solver to use
         * @return execution time of solver
         */
        double runGeneratedTest(int h, Points2D& input, Solver2D& solver);

        /**
         * Runs single test with specified input.
         *
         * @param h expected number of hull points
         * @param solver solver to use
         * @param 
         */
        double runSpecifiedTest(unsigned h, Solver2D& solver,
                                Points2D& output);

        /**
         * Asseses performance of an aprroximation scheme.
         *
         * Performance is based on relative number of hull points, since this
         * is the only parameter that concerns us when selecting solver based
         * on approximation.
         */
        void approxTests(Approximator2D& scheme);

        void runApproxInstance(Instance& inst, Approximator2D& scheme);

        /**
         * Prints single list of points
         *
         * @param points list of points to print
         * @param output target stream
         */
        void printPoints(Points2D& points, std::ofstream& output);

        void d3tests();

        /** File for logging results */
        std::ofstream logFile_;


        Points2D * globIn_;
        int counter_;

};

}
