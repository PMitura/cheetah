#pragma once

#include <cmath>

namespace ch
{

const double EPS = 1e-12;
const double PI  = 3.141592653589793238462643383279;

/** computes angle of line segment ab and x axis */
double polarAngle(double ax, double ay, double bx, double by);

/** computes distance of points a and b */
double dist(double ax, double ay, double bx, double by);

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
