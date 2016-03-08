#include "solver_2d.h"

namespace ch {

Points2D Solver2D::solve(const PointsND& inputSet)
{
    if (inputSet.getDimension() != 2)
        return Points2D();
    return algorithm(inputSet.getData());
}

}
