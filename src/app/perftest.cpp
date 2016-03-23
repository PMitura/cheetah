#include "app/perftest.h"

namespace ch
{

void PerfTest::runAllTests()
{
    // initialize tested solvers
    std::vector<Solver2D*> solvers;
    solvers.push_back(new JarvisScan2D());
    // solvers.push_back(new GrahamScan2D());
    solvers.push_back(new Quickhull2D());

    // setup parallelism
    omp_set_num_threads(1);

    // run tests
    // smallTests(solvers);

    delete *solvers.begin();
    solvers.erase(solvers.begin()); // jarvis too slow
    bigTests(solvers);

    // cleanup
    for (auto solver : solvers) {
        delete solver;
    }
}

void PerfTest::smallTests(std::vector<Solver2D*> solvers)
{
    std::vector<Instance> instances;
    instances.push_back({100 , 3,  25000, 1000});
    instances.push_back({100 , 10, 25000, 1000});
    instances.push_back({100 , 50, 25000, 1000});
    instances.push_back({1000, 3,  10000, 1000});
    instances.push_back({1000, 10, 10000, 1000});

    for (auto& inst : instances) {
        runTestInstance(inst, solvers);
    }
}

void PerfTest::bigTests(std::vector<Solver2D*> solvers)
{
    // test block instance
    std::vector<Instance> instances;
    instances.push_back({5000000, 3,     1, 1000});
    instances.push_back({5000000, 10,    1, 1000});
    instances.push_back({5000000, 50,    1, 1000});
    instances.push_back({5000000, 100,   1, 1000});
    instances.push_back({5000000, 500,   1, 1000});
    instances.push_back({5000000, 1000,  1, 1000});
    instances.push_back({5000000, 5000,  1, 1000});
    instances.push_back({5000000, 10000, 1, 1000});

    for (auto& inst : instances) {
        runTestInstance(inst, solvers);
    }
}

void PerfTest::runTestInstance(Instance& inst, std::vector<Solver2D*> solvers)
{
    std::cout << std::fixed << std::setprecision(6);
    std::cout << inst.n    << " pts, "
              << inst.h    << " on hull, "
              << inst.runs << " runs" << std::endl;
    for (auto solver : solvers) {
        std::cout << std::setw(15) << solver -> getName() << ": ";
        double timeSum = 0.0, currTime;
        int failed = 0;
        for (int i = 0; i < inst.runs; i++) {
            currTime = runGeneratedTest(inst.n, inst.h,
                    inst.span, *solver);
            if (currTime < -EPS) {
                failed++;
            } else {
                timeSum += currTime;
            }
        }
        if (failed) {
            std::cout << "[FAILED] (" << failed << " of "
                      << inst.runs << ")" << std::endl;
        } else {
            std::cout << timeSum << " ms" << std::endl;
        }
    }
    std::cout << std::endl;

}

double PerfTest::runGeneratedTest(int n, int h, double span, Solver2D& solver)
{
    Generator2D generator;
    Points2D input, output;
    generator.genUniformCircle(n, h, span, input);

    double timeStart = omp_get_wtime();
    solver.solve(input, output);
    double timeEnd   = omp_get_wtime();

    if (output.getSize() != h) {
        return -1;
    }
    return timeEnd - timeStart;
}

}
