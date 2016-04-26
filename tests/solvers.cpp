#include <gtest/gtest.h>

#include <iostream>
#include <iomanip>
#include <fstream>

#include "lib/structures.h"
#include "lib/generator.h"
#include "solvers/jarvis_scan_2d.h"
#include "solvers/graham_scan_2d.h"
#include "solvers/monotone_chain_2d.h"
#include "solvers/quickhull_2d.h"
#include "solvers/chan_2d.h"
#include "solvers/jarvis_scan_3d.h"

void printHull2D(ch::Points2D& points, std::ostream& out)
{
    const ch::data_t& data = points.getData();
    out << std::setprecision(15);
    for (auto point : data) {
        for (auto coord : point) {
            out << coord << " ";
        }
        out << std::endl;
    }
}

bool readFile(const std::string& filename, ch::Points2D& points)
{
    std::ifstream file(filename);
    if (!file.is_open())
        return false;
    long long int n;
    file >> n;
    points.clear();
    for (int i = 0; i < n; i++) {
        double a, b;
        file >> a >> b;
        points.add({a, b});
    }
    return true;
}

void testSingleFile(const std::string& filename,
                    ch::Solver2D& solver,
                    int expected)
{
    ch::Points2D input, output;
    ASSERT_TRUE(readFile(filename, input));
    solver.solve(input, output);
    // printHull2D(output, std::cout);
    EXPECT_EQ(expected, output.getSize());
}

void testSolverPremade2D(ch::Solver2D& solver)
{
    // 8 integer points, 5 on hull
    testSingleFile("tests/files/basic1.in", solver, 5);

    // 500 random points, real coords <-100, 100>, 10 on hull
    testSingleFile("tests/files/basic2.in", solver, 10);

    // 1000 random points, real coords <-0.5, 0.5>, 12 on hull
    testSingleFile("tests/files/basic3.in", solver, 12);

    // 10000 random points, real coords <-1000, 1000>, 18 on hull
    testSingleFile("tests/files/basic4.in", solver, 18);

    // 6 integer points on line
    testSingleFile("tests/files/line.in", solver, 2);

    // single point
    testSingleFile("tests/files/single.in", solver, 1);

    // diamond
    testSingleFile("tests/files/diamond.in", solver, 4);

    // hexagon
    testSingleFile("tests/files/hexagon.in", solver, 6);
}

void testSingleGen(long long n, long long h, double radius,
                   ch::Solver2D& solver)
{
    ch::Generator2D generator;
    ch::Points2D genSet, output;

    generator.genUniformCircle(n, h, radius, genSet);
    solver.solve(genSet, output);
    EXPECT_EQ(h, output.getSize());
    genSet.clear(); output.clear();

    generator.genRandomCircle(n, h, radius, genSet);
    solver.solve(genSet, output);
    EXPECT_EQ(h, output.getSize());
}

void testSolverGen2D(ch::Solver2D& solver)
{
    testSingleGen(100, 6, 1000, solver);
    testSingleGen(1000, 100, 1000, solver);
    testSingleGen(500, 10, 1, solver);
    testSingleGen(500, 250, 100000, solver);
    testSingleGen(50, 49, 100000, solver);
    testSingleGen(50, 50, 100000, solver);
}

TEST(JarvisScan2DTest, Premade)
{
    ch::JarvisScan2D jarvis;
    testSolverPremade2D(jarvis);
}

TEST(JarvisScan2DTest, Generated)
{
    ch::JarvisScan2D jarvis;
    testSolverGen2D(jarvis);
}

TEST(GrahamScan2DTest, Premade)
{
    ch::GrahamScan2D graham;
    testSolverPremade2D(graham);
}

TEST(GrahamScan2DTest, Generated)
{
    ch::GrahamScan2D graham;
    testSolverGen2D(graham);
}

TEST(MonotoneChain2DTest, Premade)
{
    ch::MonotoneChain2D mono;
    testSolverPremade2D(mono);
}

TEST(MonotoneChain2DTest, Generated)
{
    ch::MonotoneChain2D mono;
    testSolverGen2D(mono);
}

TEST(QuickHull2DTest, Premade)
{
    ch::Quickhull2D quickhull;
    testSolverPremade2D(quickhull);
}

TEST(QuickHull2DTest, Generated)
{
    ch::Quickhull2D quickhull;
    testSolverGen2D(quickhull);
}

TEST(Chan2DTest, Premade)
{
    ch::Chan2D chan;
    testSolverPremade2D(chan);
}

TEST(Chan2DTest, Generated)
{
    ch::Chan2D chan;
    testSolverGen2D(chan);
}
/*
TEST(PrintHull, EraseMe)
{
    ch::Generator2D generator;
    ch::Points2D points;

    EXPECT_TRUE(generator.genUniformCircle(1000000, 100000, 100000, points));
    std::ofstream of("example.out");
    printHull2D(points, of);
}
*/

bool readFile3D(const std::string& filename, ch::Points3D& points)
{
    std::ifstream file(filename);
    if (!file.is_open())
        return false;
    long long int n;
    file >> n;
    points.clear();
    for (int i = 0; i < n; i++) {
        double a, b, c;
        file >> a >> b >> c;
        points.add({a, b, c});
    }
    return true;
}

void testSingleFile3D(const std::string& filename, ch::Solver3D& solver)
{
    ch::Points3D input;
    ch::Polyhedron output;
    ASSERT_TRUE(readFile3D(filename, input));
    solver.solve(input, output);
}

void testSolverPremade3D(ch::Solver3D& solver)
{
    // regular tetrahedron, four triangle faces, one point inside
    testSingleFile3D("tests/files/3d_tetrahedron.in", solver);

    // regular cube, six square faces, regular tetrahedron inside
    testSingleFile3D("tests/files/3d_cube.in", solver);
}

TEST(Jarvis3DTest, Premade)
{
    ch::JarvisScan3D jarvis;
    testSolverPremade3D(jarvis);
}
