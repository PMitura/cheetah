#pragma once

#include <cmath>
#include <utility>

#include "lib/structures.h"

namespace ch
{

const double EPS = 1e-12;
const double PI  = 3.141592653589793238462643383279;

/** computes angle of line segment ab and x axis */
inline double polarAngle(double ax, double ay, double bx, double by)
{
    double res = atan2(by - ay, bx - ax);
    if (res > -EPS) {
        return res;
    } else {
        return 2*PI + res;
    }
}

/** computes dot product of two vectors */
inline double dot(double ax, double ay, double bx, double by)
{
    return ax * bx + ay * by;
}

/** computes distance of points a and b */
inline double dist(double ax, double ay, double bx, double by)
{
    return hypot(ax - bx, ay - by);
}
inline double dist(point2d_t a, point2d_t b)
{
    return hypot(a.first - b.first, a.second - b.second);
}

/** area of square given by side */
inline double square(double ax, double ay)
{
    return ax * ax + ay * ay;
}

/** difference of two vectors */
inline point2d_t vectorize(point2d_t a, point2d_t b)
{
    return {b.first - a.first, b.second - a.second};
}

/** translate point p according to vector v */
inline point2d_t translate(point2d_t p, point2d_t v)
{
    return {p.first + v.first, p.second + v.second};
}

/** scale vector by factor */
inline point2d_t scale(double vx, double vy, double factor)
{
    return {vx * factor, vy * factor};
}
inline point2d_t scale(point2d_t v, double factor)
{
    return {v.first * factor, v.second * factor};
}

/** computes line-point distance */
inline double distToLine(point2d_t a, point2d_t b, point2d_t p)
{
    point2d_t ap = vectorize(a, p),
              ab = vectorize(a, b);
    double u = dot(ap.first, ap.second, ab.first, ab.second)
               / square(ab.first, ab.second);
    point2d_t c = translate(a, scale(ab, u));
    return dist(p, c);
}

/** finds out on which side does the point lie */
int orientation(double ax, double ay,
                double bx, double by,
                double cx, double cy);

/** finds out if point is in triangle */
bool ptInTriangle(double ax, double ay,
                  double bx, double by,
                  double cx, double cy,
                  double px, double py);

}
