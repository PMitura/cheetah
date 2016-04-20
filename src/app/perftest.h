#pragma once

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <ctime>
#include <omp.h>

#include "approximators/approximator2d.h"
#include "approximators/bfp2d.h"
#include "lib/generator.h"
#include "lib/structures.h"
#include "solvers/solver_2d.h"
#include "solvers/jarvis_scan_2d.h"
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

    private:
        /** One instance of test run */
        struct Instance {
            int n, h, runs;
            double span;
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
         * Finds point where QH beats JS, and GS beats QH.
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
         * Generates, runs, and measures of single problem instance.
         *
         * Time to generate data is not counted.
         *
         * @param n size of data set
         * @param h target number of points on hull
         * @param span contraints on input set coordinates
         * @param solver solver to use
         * @return execution time of solver
         */
        double runGeneratedTest(unsigned n, unsigned h, double span, 
                                Solver2D& solver);

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

        /** File for logging results */
        std::ofstream logFile_;

};

}
