#include "lib/geometry.h"

namespace ch
{

int orientation(double ax, double ay,
                double bx, double by,
                double cx, double cy)
{
    // double cross = (ax - bx) * (ay - bx) - (ay - by) * (cy - by);
    double cross = (ax - bx) * (ay - cy) - (ay - by) * (ax - cx);

    if (fabs(cross) < EPS) return 0;  // colinear
    return (cross > EPS) ? 1 : 2; // clockwise or counter clockwise
}

bool ptInTriangle(double ax, double ay,
                  double bx, double by,
                  double cx, double cy,
                  double px, double py)
{
    int x = orientation(px, py, ax, ay, bx, by),
        y = orientation(px, py, bx, by, cx, cy),
        z = orientation(px, py, cx, cy, ax, ay);

    // exclude collinears
    if (x == 0 || y == 0 || z == 0)
        return 0;
    return ((x == y) && (y == z));
}

}
