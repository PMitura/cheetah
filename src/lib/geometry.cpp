#include "lib/geometry.h"

namespace ch
{

double polarAngle(double ax, double ay, double bx, double by)
{
    double res = atan2(by - ay, bx - ax);
    if (res > -EPS) {
        return res;
    } else {
        return 2*PI + res;
    }
}

double dist(double ax, double ay, double bx, double by)
{
    return hypot(ax - bx, ay - by);
}

int orientation(double ax, double ay,
                double bx, double by,
                double cx, double cy)
{
    double val = (ax - bx) * (cx - bx) - (ay - by) * (cy - by);
 
    if (fabs(val) < EPS) return 0;  // colinear
    return (val > EPS)? 1: 2; // clockwise or counter clockwise
}

bool ptInTriangle(double ax, double ay,
                  double bx, double by,
                  double cx, double cy,
                  double px, double py)
{
    int x = orientation(px, py, ax, ay, bx, by),
        y = orientation(px, py, bx, by, cx, cy),
        z = orientation(px, py, cx, cy, ax, ay);

    return ((x == y) && (y == z)); // exclude collinears
}

}
