#include "lib/geometry.h"

namespace ch
{

inline double polarAngle(double ax, double ay, double bx, double by)
{
    return atan2(bx - ax, by - ay);
}

}

