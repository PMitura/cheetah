#include "lib/geometry.h"

namespace ch
{

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
