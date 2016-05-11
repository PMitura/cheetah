#include "app/perftest.h"

namespace ch
{

void PerfTest::runAllTests()
{
    // initialize tested solvers
    std::vector<Solver2D*> solvers;
    // solvers.push_back(new MonotoneChain2D());
    
    // solvers.push_back(new JarvisScan2D(JarvisScan2D::POLAR));
    // solvers.push_back(new JarvisScan2D(JarvisScan2D::CROSS));

    // solvers.push_back(new GrahamScan2D());
    
    // solvers.push_back(new Quickhull2D(Quickhull2D::NAIVE));
    // solvers.push_back(new Quickhull2D(Quickhull2D::PRECOMP));
    // solvers.push_back(new Quickhull2D(Quickhull2D::FORWARD));
    
    // solvers.push_back(new Chan2D(Chan2D::GRAHAM));
    // solvers.push_back(new Chan2D(Chan2D::QUICK));
    // solvers.push_back(new Chan2D(Chan2D::COMBO));
    
    // solvers.push_back(new GrahamScan2D(GrahamScan2D::SEQ));
    // solvers.push_back(new GrahamScan2D(GrahamScan2D::PARA));
    // solvers.push_back(new GrahamScan2D(GrahamScan2D::PARA_LIN));
    // solvers.push_back(new GrahamScan2D(GrahamScan2D::PARA_STABLE));
    // solvers.push_back(new GrahamScan2D(GrahamScan2D::PARA_LIN_STABLE));
    
    solvers.push_back(new JarvisScan2D(JarvisScan2D::CROSS));
    solvers.push_back(new JarvisScan2D(JarvisScan2D::PARA));
    solvers.push_back(new JarvisScan2D(JarvisScan2D::PARA_INT));
    solvers.push_back(new JarvisScan2D(JarvisScan2D::PARA_DOUBLE));
    
    /*
    solvers.push_back(new Quickhull2D(Quickhull2D::FORWARD));
    for (int i = 1; i < 1000000; i *= 10) {
        solvers.push_back(new Quickhull2D(Quickhull2D::PARA, i));
    }
    */

    // run tests
    // smallTests(solvers);

    // delete *solvers.begin();
    // solvers.erase(solvers.begin()); // jarvis too slow
    bigTests(solvers);

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
    instances.push_back({100 , 3,  25000, 1000, 1});
    instances.push_back({100 , 10, 25000, 1000, 1});
    instances.push_back({100 , 50, 25000, 1000, 1});
    instances.push_back({1000, 3,  10000, 1000, 1});
    instances.push_back({1000, 10, 10000, 1000, 1});

    for (auto& inst : instances) {
        runTestInstance(inst, solvers);
    }
}

void PerfTest::bigTests(std::vector<Solver2D*> solvers)
{
    // test block instance
    std::vector<Instance> instances;

    /*
    for (int i = 0; i <= 1000; i += 50) {
        instances.push_back({i, std::min(i, 50), 50000, 1000, 1});
    }
    */

    instances.push_back({1000000, 100, 1, 1000, 1});
    // instances.push_back({1000000, 100, 1, 1000, 2});
    // instances.push_back({1000000, 100, 1, 1000, 4});
    // instances.push_back({1000000, 100, 1, 1000, 6});
    // instances.push_back({1000000, 100, 1, 1000, 8});
    // instances.push_back({1000000, 100, 1, 1000, 12});
    // instances.push_back({1000000, 100, 1, 1000, 16});
    // instances.push_back({1000000, 100, 1, 1000, 24});

    /*
    instances.push_back({10000000, 3,     1, 1000, 1});
    instances.push_back({10000000, 10,    1, 1000, 1});
    instances.push_back({10000000, 25,    1, 1000, 1});
    instances.push_back({10000000, 50,    1, 1000, 1});
    instances.push_back({10000000, 100,   1, 1000, 1});
    instances.push_back({10000000, 250,   1, 1000, 1});
    instances.push_back({10000000, 500,   1, 1000, 1});
    instances.push_back({10000000, 1000,  1, 1000, 1});
    instances.push_back({10000000, 5000,  1, 1000, 1});
    instances.push_back({10000000, 10000, 1, 1000, 1});
    instances.push_back({10000000, 20000, 1, 1000, 1});
    */

    globIn_ = new Points2D();
    Instance ix = instances[0];
    Generator2D generator;
    generator.genUniformCircle(ix.n, ix.h, ix.span, *globIn_);
    for (auto& inst : instances) {
        runTestInstance(inst, solvers);
    }
    delete globIn_;
}

void PerfTest::runTestInstance(Instance& inst, std::vector<Solver2D*> solvers)
{
    std::cout << std::fixed << std::setprecision(6);
    std::cout << inst.n    << " pts, "
              << inst.h    << " on hull, "
              << inst.runs << " runs" << std::endl;

    /*
    Generator2D generator;
    Points2D input;
    generator.genUniformCircle(inst.n, inst.h, inst.span, input);
    */
    omp_set_num_threads(inst.cores);
    for (auto solver : solvers) {
        // std::cout << std::setw(15) << solver -> getName() << ": ";
        double timeSum = 0.0, currTime;
        int failed = 0;
        if (solver -> getName() == "Jarvis Scan" && inst.h > 250) {
            std::cout << 0 << " s" << std::endl;
            continue;
        }
        for (int i = 0; i < inst.runs; i++) {
            currTime = runGeneratedTest(inst.h, *globIn_, *solver);
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
            std::cout << timeSum << " ";
        }
    }

    std::cout << std::endl;
}

double PerfTest::runGeneratedTest(int h, Points2D& input, Solver2D& solver)
{
    Points2D output;

    // output generated test
    // std::ofstream of("aside/compare2.in");
    // printPoints(input, of);
    // of.close();

    double timeStart = omp_get_wtime();
    solver.solve(input, output);
    double timeEnd   = omp_get_wtime();

    if ((int) output.getSize() != h) {
        return -1;
    }
    return timeEnd - timeStart;
}

void PerfTest::approxTests(Approximator2D& scheme)
{
    std::vector<Instance> instances;
    instances.push_back({5000000, 500, 1, 1000, 1});

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
        generator.genUniformCircle(inst.n, inst.h, inst.span, input);
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

void PerfTest::printPoints(Points2D& points, std::ofstream& output)
{
    output << points.getSize() << std::endl;
    output << std::fixed << std::setprecision(20);
    for (auto& pt : points.getData()) {
        output << pt[0] << " " << pt[1] << std::endl;
    }
}

}
