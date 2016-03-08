#include <gtest/gtest.h>
#include <vector>

#include "lib/cheetah.h"

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
