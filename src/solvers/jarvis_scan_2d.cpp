#include "jarvis_scan_2d.h"

namespace ch {


inline double JarvisScan2D::polarAngle(double ax, double ay, 
                                       double bx, double by) const
{
    return atan2(bx - ax, by - ay);
    /*
    if (fabs(ang) < EPS)
        return dist(pivot, a) < dist(pivot, b);
    return ang < EPS;
    */
}

Points2D& JarvisScan2D::solve(const Points2D& input, Points2D& output)
{
    std::vector<double> currPoint;
    const data_t& inputData = input.getData();

    // find point with min Y
    int minIndex = 0,
        currIndex,
        nextIndex;
    for (int i = 1; i < inputData.size(); i++) {
        if (inputData[i][1] < inputData[minIndex][1]) {
            minIndex = i;
        }
    }

    // find the rest of points
    currIndex = minIndex;
    do {
        currPoint = inputData[currIndex];
        output.add(currPoint);
        nextIndex = 0;
        for (int i = 1; i < inputData.size(); i++) {
            if (i != currIndex) {
                nextIndex = i;
            }
        }
        currIndex = nextIndex;
    } while (currIndex != minIndex);

    return output;
}

}
