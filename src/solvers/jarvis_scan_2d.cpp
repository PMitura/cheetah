#include "jarvis_scan_2d.h"

namespace ch
{

JarvisScan2D::JarvisScan2D()
{
    name_ = "Jarvis Scan";
}

Points2D& JarvisScan2D::solve(const Points2D& input, Points2D& output)
{
    // return solvePolar(input, output);
    return solveCross(input, output); // faster!
}

Points2D& JarvisScan2D::solveCross(const Points2D& input, Points2D& output)
{
    point_t currPoint;
    const data_t& inputData = input.getData();

    if (inputData.size() <= 1) {
        output = input;
        return output;
    }

    unsigned maxIndex = 0, currIndex, nextIndex;

    // find point with max Y (min X in case of tie)
    for (unsigned i = 1; i < inputData.size(); i++) {
        double dif = inputData[i][1] - inputData[maxIndex][1];
        if (fabs(dif) <= EPS) {
            if (inputData[i][0] > inputData[maxIndex][0]) {
                maxIndex = i;
            }
        } else if (dif > EPS) {
            maxIndex = i;
        }
    }

    // find the rest of points
    currIndex = maxIndex;
    do {
        output.add(inputData[currIndex]);
        // avoid setting same point as next
        nextIndex = !currIndex;

        // check orientation for all remaining n - 1 points
        for (unsigned i = 0; i < inputData.size(); i++) {
            if (i == currIndex) {
                continue;
            }

            int o = orientation(inputData[currIndex][0],
                                inputData[currIndex][1],
                                inputData[nextIndex][0],
                                inputData[nextIndex][1],
                                inputData[i        ][0],
                                inputData[i        ][1]);

            if (o == 0) {
                // exclude collinear points
                if (dist(inputData[currIndex][0],
                         inputData[currIndex][1],
                         inputData[i][0],
                         inputData[i][1])
                    >
                    dist(inputData[currIndex][0],
                         inputData[currIndex][1],
                         inputData[nextIndex][0],
                         inputData[nextIndex][1])) {
                    nextIndex = i;
                }
            } else if (o == 1) {
                // point to the left of current hull face
                nextIndex = i;
            }
        }
        currIndex = nextIndex;
    } while (currIndex != maxIndex);

    return output;
}

Points2D& JarvisScan2D::solvePolar(const Points2D& input, Points2D& output)
{
    std::vector<double> currPoint;
    const data_t& inputData = input.getData();

    if (inputData.size() <= 1) {
        output = input;
        return output;
    }

    unsigned maxIndex = 0,
             currIndex,
             nextIndex;

    // find point with max Y (min X in case of tie)
    for (unsigned i = 1; i < inputData.size(); i++) {
        double dif = inputData[i][1] - inputData[maxIndex][1];
        if (fabs(dif) <= EPS) {
            if (inputData[i][0] > inputData[maxIndex][0]) {
                maxIndex = i;
            }
        } else if (dif > EPS) {
            maxIndex = i;
        }
    }

    // find the rest of points
    currIndex = maxIndex;
    double currAngle = 0, nextAngle, minAngle, relAngle, pureAngle, nextPureAngle;
    do {
        currPoint = inputData[currIndex];
        output.add(currPoint);

        // avoid setting same point as next
        nextIndex = !currIndex;
        nextPureAngle = polarAngle(inputData[currIndex][0],
                                   inputData[currIndex][1],
                                   inputData[nextIndex][0],
                                   inputData[nextIndex][1]);
        minAngle = nextPureAngle + 2*PI - currAngle;
        if (minAngle > 2*PI + EPS) {
            minAngle -= 2*PI;
        }

        // check all remaining n - 1 points, find min polar angle
        for (unsigned i = 0; i < inputData.size(); i++) {
            if (i == currIndex) {
                continue;
            }
            pureAngle = polarAngle(inputData[currIndex][0],
                                   inputData[currIndex][1],
                                   inputData[i][0],
                                   inputData[i][1]);
            nextAngle = pureAngle + 2*PI - currAngle;
            if (nextAngle > 2*PI + EPS) {
                nextAngle -= 2*PI;
            }
            // std::cout << i << ": " << nextAngle << std::endl;
            relAngle = minAngle - nextAngle;
            if (fabs(relAngle) <= EPS) {
                // exclude collinear points
                if (dist(inputData[currIndex][0],
                         inputData[currIndex][1],
                         inputData[i][0],
                         inputData[i][1])
                    >
                    dist(inputData[currIndex][0],
                         inputData[currIndex][1],
                         inputData[nextIndex][0],
                         inputData[nextIndex][1])) {
                    minAngle = nextAngle;
                    nextPureAngle = pureAngle;
                    nextIndex = i;
                }
            } else if (relAngle > EPS) {
                minAngle = nextAngle;
                nextPureAngle = pureAngle;
                nextIndex = i;
            }
        }
        currAngle = nextPureAngle;
        currIndex = nextIndex;
    } while (currIndex != maxIndex);

    return output;
}

}
