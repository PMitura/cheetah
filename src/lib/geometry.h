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
inline double dot(point2d_t a, point2d_t b)
{
    return a.first * b.first + a.second * b.second;
}

/** computes distance of points a and b */
inline double dist(double ax, double ay, double bx, double by)
{
    return hypot(ax - bx, ay - by);
}
inline double dist(point_t& a, const point_t& b)
{
    return hypot(a[0] - b[0], a[1] - b[1]);
}
inline double dist(point2d_t a, point2d_t b)
{
    return hypot(a.first - b.first, a.second - b.second);
}

/** area of square given by side */
inline double square(point2d_t p)
{
    return p.first * p.first + p.second * p.second;
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
    double u = dot(ap, ab) / square(ab);
    point2d_t c = translate(a, scale(ab, u));
    return dist(p, c);
}

inline double partCross(const double& cx, const double& cy,
                        const double& alpha, const double& beta,
                        const double& gamma)
{
    return alpha * cx - beta * cy + gamma;
}

/** returns cross product of two vectors */
inline double cross(const double& ax, const double& ay,
                    const double& bx, const double& by,
                    const double& cx, const double& cy)
{
    return (ax - bx) * (by - cy) - (ay - by) * (bx - cx);
}

/** finds out on which side does the point lie */
inline int orientation(const double& ax, const double& ay,
                       const double& bx, const double& by,
                       const double& cx, const double& cy)
{
    double cross = (ax - bx) * (by - cy) - (ay - by) * (bx - cx);

    if (fabs(cross) < EPS) return 0;  // collinear
    return (cross > EPS) ? 1 : 2; // counter-clockwise or clockwise
}

/** finds out on which side does the point lie */
inline bool ccw(const double& ax, const double& ay,
                const double& bx, const double& by,
                const double& cx, const double& cy)
{
    return ((ax - bx) * (by - cy) - (ay - by) * (bx - cx)) > EPS;
}
inline bool ccw(const point_t& a,
                const point_t& b,
                const point_t& c)
{
    return  ((a[0] - b[0]) * (b[1] - c[1])
           - (a[1] - b[1]) * (b[0] - c[0])) > EPS;
}

/** finds out if point is in triangle */
bool ptInTriangle(double ax, double ay,
                  double bx, double by,
                  double cx, double cy,
                  double px, double py);

/** finds out if point is in given ccw ordered polygon */
inline bool ptInPolygon(const data_t& poly, const point_t& point)
{
    unsigned s = poly.size();
    if (s <= 2) {
        return 0;
    }
    for (int i = 0; i < poly.size() - 2; i++) {
        if (!ccw(poly[i], poly[i+1], poly[i+2])) {
            return 0;
        }
    }
    if (!ccw(poly[s-2], poly[s-1], poly[0])) {
        return 0;
    }
    if (!ccw(poly[s-1], poly[0], poly[1])) {
        return 0;
    }
    return 1;
}

}

