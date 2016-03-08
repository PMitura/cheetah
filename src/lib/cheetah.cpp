#include "lib/cheetah.h"

namespace ch
{

PointsND::PointsND(unsigned int dim)
    : dimension_(dim)
{}

bool PointsND::add(std::vector<double> point)
{
    if (point.size() != dimension_)
        return false;
    data_.push_back(point);
    return true;
}

Points2D::Points2D()
    : PointsND(2)
{}

}
