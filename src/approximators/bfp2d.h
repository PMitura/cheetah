#pragma once

#include <algorithm>

#include "approximators/approximator2d.h"
#include "lib/geometry.h"
#include "lib/structures.h"

namespace ch
{

class BFP2D : public Approximator2D
{
    public:
        BFP2D();
        Points2D& approximate(const Points2D& input, Points2D& output);
        int maxReachable() const;

    private:
        int stripsCount_;
        Points2D& sequential(const Points2D& input, Points2D& output);

};

}
