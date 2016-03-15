#include "solvers/quickhull_2d.h"

namespace ch
{

Quickhull2D::Quickhull2D()
    :globOut_(NULL)
{}

Points2D& Quickhull2D::solve(const Points2D& input, Points2D& output)
{
    return solveNaive(input, output);
}

void Quickhull2D::recNaive(point_t& a, point_t& b, data_t* plane)
{
    if ((*plane).size() == 0)
        return;

    // find point furthest from ab
    for (auto& pt : (*plane)) {

    }
}

Points2D& Quickhull2D::solveNaive(const Points2D& input, Points2D& output)
{
    // find min and max X
    
}

}
