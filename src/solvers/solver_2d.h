#pragma once

#include <vector>
#include <lib/cheetah.h>

namespace ch
{

typedef std::vector<std::pair<int, int>> points2D;

class Solver2D
{
    public:
        virtual Points2D solve(const Points2D& inputSet) = 0;
};

}
