#include <gtest/gtest.h>

#include "lib/geometry.h"

TEST(GeometryTest, PointToLine)
{
    EXPECT_DOUBLE_EQ(5.0, ch::distToLine({0, 0}, {0, 5}, {5, 5}));
    EXPECT_DOUBLE_EQ(sqrt(0.5), ch::distToLine({0, 0}, {1, 1}, {0, 1}));
}
