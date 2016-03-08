#include <gtest/gtest.h>

#include "lib/cheetah.h"

TEST(PointsNDTest, Constructor) {
    ch::PointsND a(66);
    EXPECT_EQ(66, a.getDimension());
}
