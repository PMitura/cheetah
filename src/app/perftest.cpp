#include "app/perftest.h"

namespace ch
{

void PerfTest::runAllTests()
{
    // initialize tested solvers
    std::vector<Solver2D*> solvers;
    solvers.push_back(new JarvisScan2D());
    solvers.push_back(new GrahamScan2D());
    solvers.push_back(new Quickhull2D());

    // cleanup
    for (auto& solver : solvers) {
        delete solver;
    }
}

bool PerfTest::runGeneratedTest(int n, int h, double span, Solver2D& solver)
{
    Generator2D generator;
    Points2D input, output;
    generator.genRandomCircle(n, h, span, input);

    double timeStart = omp_get_wtime();
    solver.solve(input, output);
    double timeEnd   = omp_get_wtime();

    if (output.getSize() != h) {
        std::cout << "FAILED" << std::endl;
        return false;
    }
    std::cout << timeEnd - timeStart << std::endl;
    return true;
}

}
