#include "jarvis_scan_2d.h"

namespace ch {

Points2D& JarvisScan2D::solve(const Points2D& input, Points2D& output)
{
    std::vector<double> currPoint;
    const data_t& inputData = input.getData();

    // find point with min Y
    unsigned minIndex = 0,
             currIndex,
             nextIndex;
    for (unsigned i = 1; i < inputData.size(); i++) {
        if (inputData[i][1] < inputData[minIndex][1]) {
            minIndex = i;
        }
    }

    // find the rest of points
    currIndex = minIndex;
    do {
        currPoint = inputData[currIndex];
        output.add(currPoint);

        // avoid checking with self
        nextIndex = !currIndex;
        double maxAngle = polarAngle(inputData[currIndex][0],
                                     inputData[currIndex][1],
                                     inputData[nextIndex][0],
                                     inputData[nextIndex][1]);

        // check all n - 1 points, find max polar angle
        double nextAngle, relAngle;
        for (unsigned i = nextIndex + 1; i < inputData.size(); i++) {
            if (i == currIndex) {
                continue;
            }
            nextAngle = polarAngle(inputData[currIndex][0],
                                   inputData[currIndex][1],
                                   inputData[i][0],
                                   inputData[i][1]);
            relAngle = nextAngle - maxAngle;
            if (fabs(relAngle) < EPS) {
                // collinear points
                if (dist(inputData[currIndex][0],
                         inputData[currIndex][1],
                         inputData[i][0],
                         inputData[i][1])
                    <
                    dist(inputData[currIndex][0],
                         inputData[currIndex][1],
                         inputData[nextIndex][0],
                         inputData[nextIndex][1])) {
                    maxAngle = nextAngle;
                    nextIndex = i;
                }
            } else if (relAngle < EPS) {
                maxAngle = nextAngle;
                nextIndex = i;
            }
        }
        currIndex = nextIndex;
    } while (currIndex != minIndex);

    return output;
}

}
