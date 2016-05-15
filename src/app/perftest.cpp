#include "app/perftest.h"

namespace ch
{

void PerfTest::runAllTests()
{
    // initialize tested solvers
    std::vector<Solver2D*> solvers;

    // d3tests();
    // solvers.push_back(new MonotoneChain2D());
    
    // solvers.push_back(new JarvisScan2D(JarvisScan2D::POLAR));
    solvers.push_back(new JarvisScan2D(JarvisScan2D::CROSS));
    
    // solvers.push_back(new Quickhull2D(Quickhull2D::NAIVE));
    // solvers.push_back(new Quickhull2D(Quickhull2D::PRECOMP));
    // solvers.push_back(new Quickhull2D(Quickhull2D::FORWARD));
    
    // solvers.push_back(new Chan2D(Chan2D::GRAHAM));
    // solvers.push_back(new Chan2D(Chan2D::QUICK));
    // solvers.push_back(new Chan2D(Chan2D::COMBO));
    
    // solvers.push_back(new JarvisScan2D(JarvisScan2D::CROSS));
    // solvers.push_back(new JarvisScan2D(JarvisScan2D::PARA));
    // solvers.push_back(new JarvisScan2D(JarvisScan2D::PARA_INT));
    // solvers.push_back(new JarvisScan2D(JarvisScan2D::PARA_DOUBLE));

    // solvers.push_back(new GrahamScan2D(GrahamScan2D::SEQ));
    // solvers.push_back(new GrahamScan2D(GrahamScan2D::PARA));
    // solvers.push_back(new GrahamScan2D(GrahamScan2D::PARA_LIN));
    // solvers.push_back(new GrahamScan2D(GrahamScan2D::PARA_STABLE));
    // solvers.push_back(new GrahamScan2D(GrahamScan2D::PARA_LIN_STABLE));
    
    /*
    solvers.push_back(new Quickhull2D(Quickhull2D::FORWARD));
    for (int i = 1; i < 1000000; i *= 10) {
        solvers.push_back(new Quickhull2D(Quickhull2D::PARA, i));
    }
    */
    // solvers.push_back(new Quickhull2D(Quickhull2D::PARA, 10000));
    
    // solvers.push_back(new Chan2D(Chan2D::COMBO));
    // solvers.push_back(new Chan2D(Chan2D::PARA_OVER));
    // solvers.push_back(new Chan2D(Chan2D::PARA_ALGO));
    // solvers.push_back(new Chan2D(Chan2D::PARA_COMBO));

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

    for (int i = 0; i <= 2000000; i += 250000) {
        instances.push_back({i, std::min(i, 200), 1, 1000, 1});
    }

    /*
    instances.push_back({5000000, 100, 1, 1000, 1});
    instances.push_back({5000000, 100, 1, 1000, 2});
    instances.push_back({5000000, 100, 1, 1000, 4});
    instances.push_back({5000000, 100, 1, 1000, 6});
    instances.push_back({5000000, 100, 1, 1000, 8});
    instances.push_back({5000000, 100, 1, 1000, 12});
    instances.push_back({5000000, 100, 1, 1000, 16});
    instances.push_back({5000000, 100, 1, 1000, 24});
    */

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

    /*
    globIn_ = new Points2D();
    Instance ix = instances[0];
    Generator2D generator;
    generator.genUniformCircle(ix.n, ix.h, ix.span, *globIn_);
    */
    
    counter_ = 0;
    for (auto& inst : instances) {
        runTestInstance(inst, solvers);
    }
}

double PerfTest::runTest(int n, int h, int span, int runs,
                int cores, SolverType type)
{
    Solver2D * solver = NULL;
    if (cores > 1) {
        switch (type) {
            case JARVIS:
                solver = new JarvisScan2D(JarvisScan2D::PARA_DOUBLE);
                break;
            case GRAHAM:
                solver = new GrahamScan2D(GrahamScan2D::PARA);
                break;
            case QUICKHULL:
                solver = new Quickhull2D(Quickhull2D::PARA);
                break;
            case CHAN:
                solver = new Chan2D(Chan2D::PARA_OVER);
                break;
            case ANDREW:
                solver = new MonotoneChain2D();
                break;
            default:
                solver = new Quickhull2D(Quickhull2D::PARA);
                break;
        }
    } else {
        switch (type) {
            case JARVIS:
                solver = new JarvisScan2D();
                break;
            case GRAHAM:
                solver = new GrahamScan2D();
                break;
            case QUICKHULL:
                solver = new Quickhull2D();
                break;
            case CHAN:
                solver = new Chan2D();
                break;
            case ANDREW:
                solver = new MonotoneChain2D();
                break;
            default:
                solver = new Quickhull2D();
                break;
        }
    }

    Points2D input;
    Generator2D generator;
    std::cout << "Generating test instance..." << std::endl;
    generator.genUniformCircle(n, h, span, input);
    std::cout << "...done, running test" << std::endl;
    omp_set_num_threads(cores);
    double timeSum = 0.0, currTime;
    int failed = 0;
    for (int i = 0; i < runs; i++) {
        currTime = runGeneratedTest(h, input, *solver);
        if (currTime < -EPS) {
            failed++;
        } else {
            timeSum += currTime;
        }
    }
    delete solver;
    return timeSum;
}

void PerfTest::runTestInstance(Instance& inst, std::vector<Solver2D*> solvers)
{
    std::cout << std::fixed << std::setprecision(6);
    std::cout << inst.n    << " pts, "
              << inst.h    << " on hull, "
              << inst.runs << " runs" << std::endl;

    Generator2D generator;
    Points2D input;
    generator.genUniformCircle(inst.n, inst.h, inst.span, input);

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
            currTime = runGeneratedTest(inst.h, input, *solver);
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
        std::cout << "[WARNING] Precision errors occured (h may be too high)"
            << std::endl;
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

void PerfTest::d3tests() {
    std::ifstream ifs("aside/3d.data");
    int n;
    ifs >> n;
    Points3D input;
    for (int i = 0; i < n; i++) {
        double a, b, c;
        ifs >> a >> b >> c;
        point_t pt = {a, b, c};
        input.add(pt);
    }
    JarvisScan3D jarv;
    Polyhedron output;
    double timeA = omp_get_wtime();
    jarv.solve(input, output);
    double timeB = omp_get_wtime();
    std::cout << timeB - timeA << std::endl;
    std::cout << "faces: " << output.getFaces().size() << std::endl;
}

}
