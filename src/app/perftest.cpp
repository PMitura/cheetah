#include "app/perftest.h"

namespace ch
{

void PerfTest::runAllTests()
{
    // initialize tested solvers
    std::vector<Solver2D*> solvers;
    solvers.push_back(new GrahamScan2D());
    solvers.push_back(new MonotoneChain2D());
    solvers.push_back(new Quickhull2D());
    solvers.push_back(new JarvisScan2D());

    // setup parallelism
    omp_set_num_threads(1);

    // run tests
    // smallTests(solvers);

    // delete *solvers.begin();
    // solvers.erase(solvers.begin()); // jarvis too slow
    // bigTests(solvers);
    testEdges();

    // approximation tests
    // BFP2D bfp; approxTests(bfp);

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

    /*
    instances.push_back({5000000, 3,     1, 1000});
    instances.push_back({5000000, 10,    1, 1000});
    instances.push_back({5000000, 50,    1, 1000});
    instances.push_back({5000000, 100,   1, 1000});
    instances.push_back({5000000, 500,   1, 1000});
    instances.push_back({5000000, 1000,  1, 1000});
    instances.push_back({5000000, 5000,  1, 1000});
    instances.push_back({5000000, 10000, 1, 1000});
    */

    instances.push_back({1000000, 3,  1, 1000});
    instances.push_back({1000000, 6,  1, 1000});
    instances.push_back({1000000, 10, 1, 1000});
    instances.push_back({1000000, 17, 1, 1000});
    for (int hull = 25; hull <= 1000; hull += 25) {
        instances.push_back({1000000, hull, 1, 1000});
    }

    logFile_.open("aside/log2.data");
    if (!logFile_.is_open()) {
        std::cout << "Error opening logfile" << std::endl;
        return;
    }
    for (auto& inst : instances) {
        // don't run jarvis later
        if (inst.h == 200) {
            delete *solvers.end();
            solvers.pop_back();
        }
        runTestInstance(inst, solvers);
    }
    logFile_.close();
}

void PerfTest::testEdges()
{
    logFile_.open("aside/edge2.data");
    if (!logFile_.is_open()) {
        std::cout << "Error opening logfile" << std::endl;
        return;
    }
    for (int i = 10; i < 100000; i += 10000) {
        findEdge(i);
    }
}

void PerfTest::findEdge(int setSize)
{
    JarvisScan2D jarvis;
    GrahamScan2D graham;
    Quickhull2D quick;
    logFile_ << setSize;
    std::cout << std::endl << "Find edge for n = " << setSize << std::endl;

    int hull = 3;
    double tJarvis = 0, tGraham = 1e9, tQuick = 1;
    while (tJarvis < tQuick) {
        tJarvis = runGeneratedTest(setSize, hull, 1000, jarvis);
        tQuick  = runGeneratedTest(setSize, hull, 1000, quick);
        hull++;
    }
    logFile_ << " " << hull;
    std::cout << "  Jarvis - Quick: " << hull << std::endl;

    while (tQuick < tGraham) {
        tQuick  = runGeneratedTest(setSize, hull, 1000, quick);
        tGraham = runGeneratedTest(setSize, hull, 1000, graham);
        hull += 10;
    }
    logFile_ << " " << hull << std::endl;
    std::cout << "  Quick - Graham: " << hull << std::endl;
}

void PerfTest::runTestInstance(Instance& inst, std::vector<Solver2D*> solvers)
{
    std::cout << std::fixed << std::setprecision(6);
    std::cout << inst.n    << " pts, "
              << inst.h    << " on hull, "
              << inst.runs << " runs" << std::endl;
    if (logFile_.is_open()) {
        logFile_ << inst.h;
    }

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
        if (logFile_.is_open()) {
            logFile_ << " " << timeSum;
        }
    }

    std::cout << std::endl;
    if (logFile_.is_open()) {
        logFile_ << std::endl;
    }
}

double PerfTest::runGeneratedTest(unsigned n, unsigned h, double span,
                                  Solver2D& solver)
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

void PerfTest::approxTests(Approximator2D& scheme)
{
    std::vector<Instance> instances;
    instances.push_back({5000000, 500, 1, 1000});

    for (auto& inst : instances) {
        runApproxInstance(inst, scheme);
    }
}

void PerfTest::runApproxInstance(Instance& inst, Approximator2D& scheme)
{
    std::cout << std::endl << "Approximation " << scheme.getName()
              << " on set of size " << inst.n << " points, " << inst.h
              << " on hull, averaged from " << inst.runs << " runs."
              << std::endl;
    int totalHull = 0;
    double totalTime = 0;
    Generator2D generator;
    Points2D input, output;
    for (int i = 0; i < inst.runs; i++) {
        output.clear();
        generator.genRandomCircle(inst.n, inst.h, inst.span, input);
        double timeStart = omp_get_wtime();
        scheme.approximate(input, output);
        double timeEnd = omp_get_wtime();

        totalHull += output.getSize();
        totalTime += timeEnd - timeStart;
    }
    double reached = ((double) totalHull)
                   / ((double) scheme.maxReachable() * inst.runs),
           target  = ((double) inst.h) / inst.n;

    std::cout << "  avg on hull:    " << (double) totalHull / inst.runs
              << " points" << std::endl;
    std::cout << "  target percent: " << target * 100 << "%"
              << std::endl;
    std::cout << "  actual percent: " << reached * 100 << "%"
              << std::endl;
    std::cout << "  total time " << totalTime << " ms" << std::endl;
}

}
