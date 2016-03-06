#pragma once

#include <vector>

namespace ch
{

typedef std::vector<std::pair<int, int>> points2D;

class Solver2D
{
    public:
        virtual points2D solve(points2D inputSet) = 0;
};

}
