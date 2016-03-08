#include "lib/cheetah.h"

namespace ch
{

PointsND::PointsND(int dim)
    : dimension_(dim)
{}

Points2D::Points2D()
    : PointsND(2)
{}

}
