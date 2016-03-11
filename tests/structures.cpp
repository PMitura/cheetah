#include <gtest/gtest.h>
#include <vector>

#include "lib/structures.h"
#include "lib/generator.h"

TEST(PointsNDTest, Constructor)
{
    ch::PointsND a(66), b(-1);
    EXPECT_EQ(66, a.getDimension());
    EXPECT_EQ(-1, b.getDimension());
}

TEST(PointsNDTest, PointAddition)
{
    ch::PointsND a(3);
    std::vector<double> x = {1, 2, 3},
                        y = {0.1, 0, -5.0},
                        z = {1, 2, 3, 4};
    a.add(x);
    a.add(y);
    EXPECT_EQ(2, a.getSize());
    EXPECT_FALSE(a.add(z));
    EXPECT_EQ(2, a.getSize());
}

TEST(Points2DTest, Constructor)
{
    ch::Points2D a;
    EXPECT_EQ(2, a.getDimension());
}

TEST(Generator2DTestUniform, EdgeCases)
{
    ch::Generator2D generator;
    ch::Points2D points;
    EXPECT_FALSE(generator.genUniformCircle(0, 0, 10, points));
    EXPECT_FALSE(generator.genUniformCircle(5, 10, 10, points));
    EXPECT_FALSE(generator.genUniformCircle(10, 4, 0, points));
    EXPECT_FALSE(generator.genUniformCircle(1, 1, 10, points));
    EXPECT_EQ(0, points.getSize());

    EXPECT_TRUE(generator.genUniformCircle(100, 10, 100, points));
    EXPECT_EQ(100, points.getSize());
    points.clear();
    EXPECT_TRUE(generator.genUniformCircle(500, 250, 100, points));
    EXPECT_EQ(500, points.getSize());
    points.clear();
    EXPECT_TRUE(generator.genUniformCircle(100000, 50000, 100000, points));
    EXPECT_EQ(100000, points.getSize());
}

TEST(Generator2DTestRandom, EdgeCases)
{
    ch::Generator2D generator;
    ch::Points2D points;

    EXPECT_TRUE(generator.genRandomCircle(100, 10, 100, points));
    EXPECT_EQ(100, points.getSize());
    points.clear();
    EXPECT_TRUE(generator.genRandomCircle(500, 250, 100, points));
    EXPECT_EQ(500, points.getSize());
    points.clear();
    EXPECT_TRUE(generator.genRandomCircle(100000, 50000, 100000, points));
    EXPECT_EQ(100000, points.getSize());
}

