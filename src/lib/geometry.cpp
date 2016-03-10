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

}
