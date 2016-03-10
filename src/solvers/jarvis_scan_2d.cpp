#include "jarvis_scan_2d.h"

namespace ch {

Points2D& JarvisScan2D::solve(const Points2D& input, Points2D& output)
{
    std::vector<double> currPoint;
    const data_t& inputData = input.getData();

    // find point with max Y (min X in case of tie)
    unsigned minIndex = 0,
             currIndex,
             nextIndex;
    for (unsigned i = 1; i < inputData.size(); i++) {
        double dif = inputData[i][1] - inputData[minIndex][1];
        if (fabs(dif) <= EPS) {
            if (inputData[i][0] > inputData[minIndex][0]) {
                minIndex = i;
            }
        } else if (dif > EPS) {
            minIndex = i;
        }
    }

    // find the rest of points
    currIndex = minIndex;
    double currAngle = 0, nextAngle, relAngle;
    do {
        std::cout << "pt " << inputData[currIndex][0] << ", " << inputData[currIndex][1] << std::endl;
        currPoint = inputData[currIndex];
        output.add(currPoint);

        // avoid setting same point as next
        nextIndex = !currIndex;
        double minAngle = polarAngle(inputData[currIndex][0],
                                     inputData[currIndex][1],
                                     inputData[nextIndex][0],
                                     inputData[nextIndex][1])
                          + currAngle;

        // check all n - 1 points, find min polar angle
        for (unsigned i = 0; i < inputData.size(); i++) {
            if (i == currIndex) {
                continue;
            }
            nextAngle = polarAngle(inputData[currIndex][0],
                                   inputData[currIndex][1],
                                   inputData[i][0],
                                   inputData[i][1])
                        + currAngle;
            if (nextAngle > 2*PI - EPS)
                nextAngle -= 2*PI;
            std::cout << i << ": " << nextAngle << std::endl;
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
                    nextIndex = i;
                }
            } else if (relAngle > EPS) {
                minAngle = nextAngle;
                nextIndex = i;
            }
        }
        currAngle = nextAngle;
        currIndex = nextIndex;
    } while (currIndex != minIndex);

    return output;
}

}
