#pragma once

#include "lib/cheetah.h"
#include "solvers/solver_2d.h"

namespace ch {

class JarvisScan2D : public Solver2D
{
    protected:
        Points2D algorithm(const data_t& input);
};

}
