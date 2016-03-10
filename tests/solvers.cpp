#include <gtest/gtest.h>

#include <iostream>
#include <fstream>

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

bool readFile(std::string filename, ch::Points2D& points)
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

void testSolver2D(ch::Solver2D& solver)
{
    ch::Points2D input, output;

    ASSERT_TRUE(readFile("tests/files/basic1.in", input));
    solver.solve(input, output);
    EXPECT_EQ(5, output.getSize());
    output.clear();
    ASSERT_TRUE(readFile("tests/files/basic2.in", input));
    solver.solve(input, output);
    EXPECT_EQ(10, output.getSize());
    output.clear();
    ASSERT_TRUE(readFile("tests/files/basic3.in", input));
    solver.solve(input, output);
    EXPECT_EQ(12, output.getSize());
    output.clear();
    ASSERT_TRUE(readFile("tests/files/basic4.in", input));
    solver.solve(input, output);
    EXPECT_EQ(18, output.getSize());
}

TEST(JarvisScan2DTest, PremadeData)
{
    ch::JarvisScan2D jarvis;
    testSolver2D(jarvis);
}
