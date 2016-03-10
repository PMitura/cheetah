#include <gtest/gtest.h>

#include <iostream>

#include "lib/structures.h"
#include "solvers/jarvis_scan_2d.h"

void printHull2D(ch::Points2D& points)
{
    const ch::data_t& data = points.getData();
    for (auto point : data) {
        for (auto coord : point) {
            std::cout << coord << " ";
        }
        std::cout << std::endl;
    }
}

void testSolver2D(ch::Solver2D& solver)
{
    ch::Points2D small, smallResult;
    small.add({0, 0});
    small.add({1, 0});
    small.add({0, 1});
    small.add({2, 0});
    small.add({0, 2});
    small.add({2, 2});
    small.add({1, 3});
    small.add({1, 1});
    solver.solve(small, smallResult);
    EXPECT_EQ(5, smallResult.getSize());
    printHull2D(smallResult);
}

TEST(JarvisScan2DTest, PremadeData)
{
    ch::JarvisScan2D jarvis;
    testSolver2D(jarvis);
}
