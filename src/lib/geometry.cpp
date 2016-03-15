#include "lib/geometry.h"

namespace ch
{

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
