#include <iostream>
#include <fstream>
#include <unistd.h>
#include <omp.h>

#include "app/perftest.h"

#include "lib/lib.h"
#include "lib/structures.h"

#include "solvers/chan_2d.h"
#include "solvers/graham_scan_2d.h"
#include "solvers/jarvis_scan_2d.h"
#include "solvers/jarvis_scan_3d.h"
#include "solvers/monotone_chain_2d.h"
#include "solvers/quickhull_2d.h"
#include "solvers/solver_2d.h"

using namespace std;

#ifndef TEST_MAIN

bool read2D(std::istream& is, ch::Points2D& input)
{
    long long int n = 0;
    if (!(is >> n)) {
        std::cerr << "[ERROR] Cannot read input size" << std::endl;
        return 0;
    }
    if (n > (1LL << 31)) {
        std::cerr << "[ERROR] Input is too large (max 2^31)" << std::endl;
        return 0;
    }

    double a, b;
    for (int i = 0; i < n; i++) {
        if (!(is >> a >> b)) {
            std::cerr << "[ERROR] Cannot read point" << std::endl;
            return 0;
        }
        input.add({a, b});
    }

    return 1;
}

void write2D(std::ostream& os, ch::Points2D& output)
{
    const std::vector<std::vector<double>>& d = output.getData();

    os << d.size() << std::endl;
    for (auto& i : d) {
        os << i[0] << " " << i[1] << std::endl;
    }
}

bool read3D(std::istream& is, ch::Points3D& input)
{
    long long int n;
    is >> n;
    if (n > (1LL << 31)) {
        std::cerr << "[ERROR] Input is too large (max 2^31)" << std::endl;
        return 0;
    }

    double a, b, c;
    for (int i = 0; i < n; i++) {
        if (!(is >> a >> b >> c)) {
            std::cerr << "[ERROR] Cannot read points" << std::endl;
            return 0;
        }
        input.add({a, b, c});
    }

    return 1;
}

void write3D(std::ostream& os, ch::Polyhedron& output)
{
    const std::vector<ch::Points3D>& d = output.getFaces();

    os << d.size() << std::endl;
    for (auto& i : d) {
        const std::vector<std::vector<double>>& f = i.getData();
        os << f.size() << std::endl;
        for (auto& j : f) {
            os << j[0] << " " << j[1] << " " << j[2] << std::endl;
        }
    }
}

/**
 * Interface for command line application
 *
 * Kind of spaghetti
 */
int main(int argc, char ** argv)
{
    // if custom input
    bool useFileInput = 0;
    string inputFilename;
    ifstream inputFile;

    // if custom output
    bool useFileOutput = 0;
    string outputFilename;
    ifstream outputFile;

    // perftest related
    bool perfMode = 0;
    std::vector<long long int> instance;
    int index;
    
    // if custom algo
    bool customAlgo = 0;
    string usedAlgo;

    int dimension = 2;

    bool displayTime = 0;

    bool parallelMode = 0;
    int threads = 0;

    int c, val;
    bool endFlag = 0, wasError = 0, argOK;
    string nxtArg;
    while ((c = getopt(argc, argv, "i:o:p:s:d:m:t"))) {
        switch (c) {
            case 'i':
                useFileInput = 1;
                inputFilename = optarg;
                R("got input " << inputFilename);
                break;

            case 'o':
                useFileOutput = 1;
                outputFilename = optarg;
                R("got output " << outputFilename);
                break;

            case 'p':
                index = optind - 1;
                perfMode = 1;
                while (index < argc) {
                    nxtArg = argv[index++];
                    if (nxtArg[0] == '-') {
                        break;
                    }
                    argOK = 1;
                    val = 0;
                    for (auto& c : nxtArg) {
                        if (!isdigit(c)) {
                            argOK = 0;
                            break;
                        }
                        val *= 10;
                        val += c - '0';
                    }
                    if (!argOK) {
                        std::cerr << "[ERROR] Non-numerical or negative "
                            << "argument in perftest" << std::endl;
                        endFlag = 1;
                        wasError = 1;
                        break;
                    }
                    R("instance value " << val);
                    instance.push_back(val);
                }
                if (endFlag) {
                    break;
                }
                if (instance.size() != 5) {
                    std::cerr << "[ERROR] Incorrect number of arguments "
                              << "for performance testing" << std::endl;
                    wasError = 1;
                }
                optind = index - 1;
                break;

            case 's':
                customAlgo = 1;
                usedAlgo = optarg;
                R("custom algo: " << usedAlgo);
                break;

            case 'd':
                nxtArg = optarg;
                if (nxtArg == "2") {
                    break;
                }
                if (nxtArg == "3") {
                    R("3D algo");
                    dimension = 3;
                    break;
                }
                std::cerr << "[ERROR] Wrong dimension: " << nxtArg << std::endl;
                endFlag = 1;
                wasError = 1;
                break;

            case 't':
                R("timing enabled");
                displayTime = 1;
                break;

            case 'm':
                R("set parellel mode");
                parallelMode = 1;
                nxtArg = optarg;
                for (auto& c : nxtArg) {
                    if (!isdigit(c)) {
                        argOK = 0;
                        break;
                    }
                    threads *= 10;
                    threads += c - '0';
                }
                if (!argOK) {
                    std::cerr << "[ERROR] Non-numerical or negative"
                        << " number of threads" << std::endl;
                    return 1;
                }
                if (threads < 1) {
                    std::cerr << "[ERROR] Cannot use zero threads" << std::endl;
                    return 1;
                }
                if (threads > (int) ch::MAX_NUM_THREADS) {
                    std::cerr << "[ERROR] Too many threads " <<
                        "(limit set at " << ch::MAX_NUM_THREADS << ")" 
                        << std::endl;
                    return 1;
                }

                break;

            case ':':
                std::cerr << "[ERROR] Missing option argument" << std::endl;
                endFlag = 1;
                wasError = 1;
                break;

            case '?':
                std::cerr << "[ERROR] Unknown option" << std::endl;
                endFlag = 1;
                wasError = 1;
                break;

            default:
                endFlag = 1;
        }
        if (endFlag)
            break;
    }
    if (wasError) {
        return 1;
    } 

    if ((useFileInput || useFileOutput) && perfMode) {
        std::cout << "[WARNING] Performance testing does not use any I/O"
            << std::endl;
    }

    if (dimension == 2) {
        ch::SolverType sType = ch::QUICKHULL;
        if (customAlgo) {
            if (usedAlgo == "jarvis") {
                sType = ch::JARVIS;
            } else if (usedAlgo == "graham") {
                sType = ch::GRAHAM;
            } else if (usedAlgo == "quickhull") {
                sType = ch::QUICKHULL;
            } else if (usedAlgo == "chan") {
                sType = ch::CHAN;
            } else if (usedAlgo == "andrew") {
                // experimental
                sType = ch::ANDREW;
            } else {
                std::cerr << "[ERROR] Unknown 2D solver type: " << usedAlgo <<
                    std::endl;
                return 1;
            }
        }

        if (perfMode) {
            R("Running perftest");
            if (instance[0] > (1LL << 31)) {
                std::cerr << "[ERROR] Too large instance" <<  std::endl;
                return 1;
            }
            if (instance[1] > instance[0]) {
                std::cerr << "[ERROR] Selected number of points on hull "
                    << "is greater than number of points in set" << 
                    std::endl;
                return 1;
            }
            if (instance[4] > (int) ch::MAX_NUM_THREADS) {
                std::cerr << "[ERROR] Number of threads higher than limit"
                    << std::endl;
                return 1;
            } else if (instance[4] < 1) {
                std::cerr << "[ERROR] Number of threads lesser than one"
                    << std::endl;
                return 1;
            }

            ch::PerfTest ptest;
            double perf = ptest.runTest(instance[0], instance[1],
                    instance[2], instance[3], instance[4], sType);
            std::cout << "time: " << perf << std::endl;
            return 0;
        }

        ch::Points2D input;
        if (useFileInput) {
            std::ifstream ifs(inputFilename);
            if (!ifs.is_open()) {
                std::cerr << "[ERROR] Cannot open input file " << inputFilename 
                    << std::endl;
                return 1;
            }
            if (!read2D(ifs, input)) {
                return 1;
            }
            ifs.close();
        } else {
            if (!read2D(std::cin, input)) {
                return 1;
            }
        }

        ch::Points2D output;
        double timeA = omp_get_wtime();
        if (parallelMode) {
            D("parallel " << threads);
            findHullParallel(input, output, sType, threads);
        } else {
            findHull(input, output, sType);
        }
        double timeB = omp_get_wtime();
        if (displayTime) {
            std::cout << "Execution time: " << timeB - timeA << " ms." << std::endl;
        }

        if (useFileOutput) {
            std::ofstream ofs(outputFilename);
            if (!ofs.is_open()) {
                std::cerr << "[ERROR] Cannot open output file " << outputFilename 
                    << std::endl;
                return 1;
            }
            write2D(ofs, output);
        } else {
            write2D(std::cout, output);
        }
    } else {
        // three dimensions
        if (perfMode) {
            std::cerr << "[ERROR] Performance testing not implemented for 3D: "
                << usedAlgo << std::endl;
            return 1;
        }

        if (parallelMode) {
            std::cout << "[WARNING] No 3D parallel algorithm available, "
                << "using sequential solver" << std::endl;
        }

        if (customAlgo && usedAlgo != "jarvis") {
            std::cerr << "[ERROR] Unknown 3D solver type: " << usedAlgo <<
                std::endl;
            return 1;
        }

        ch::Points3D input;
        if (useFileInput) {
            std::ifstream ifs(inputFilename);
            if (!ifs.is_open()) {
                std::cerr << "[ERROR] Cannot open input file " << inputFilename 
                    << std::endl;
                return 1;
            }
            if (!read3D(ifs, input)) {
                return 1;
            }
            ifs.close();
        } else {
            if (!read3D(std::cin, input)) {
                return 1;
            }
        }

        ch::Polyhedron output;
        double timeA = omp_get_wtime();
        findHull3D(input, output);
        double timeB = omp_get_wtime();
        if (displayTime) {
            std::cout << "Execution time: " << timeB - timeA << " ms." << std::endl;
        }

        if (useFileOutput) {
            std::ofstream ofs(outputFilename);
            if (!ofs.is_open()) {
                std::cerr << "[ERROR] Cannot open output file " << outputFilename 
                    << std::endl;
                return 1;
            }
            write3D(ofs, output);
        } else {
            write3D(std::cout, output);
        }
    }

    return 0;
}

#endif
