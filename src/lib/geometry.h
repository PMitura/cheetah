#pragma once

#include <cmath>

namespace ch
{

const double EPS = 1e-12;

/** computes angle of line segment ab and x axis */
inline double polarAngle(double ax, double ay, double bx, double by)
{
    return atan2(bx - ax, by - ay);
}

/** computes distance of points a and b */
inline double dist(double ax, double ay, double bx, double by)
{
    return hypot(ax - bx, ay - by);
}

}
