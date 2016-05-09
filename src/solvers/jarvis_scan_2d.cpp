#include "jarvis_scan_2d.h"

namespace ch
{

JarvisScan2D::JarvisScan2D()
{
    name_ = "Jarvis Scan";
    variant_ = CROSS;
}

JarvisScan2D::JarvisScan2D(Variant v)
{
    name_ = "Jarvis Scan";
    variant_ = v;
}


Points2D& JarvisScan2D::solve(const Points2D& input, Points2D& output)
{
    switch (variant_) {
        case POLAR:
            return solvePolar(input, output);
        case CROSS:
            return solveCross(input, output);
        case PARA:
            return solvePara(input, output);
    }
    // should never happen
    return solveCross(input, output);
}

Points2D& JarvisScan2D::solveCross(const Points2D& input, Points2D& output)
{
    if (input.getSize() <= 2) {
        output = input;
        return output;
    }

    point_t currPoint;
    const data_t& inputData = input.getData();

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
            } else if (o == 2) {
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
    if (input.getSize() <= 1) {
        output = input;
        return output;
    }

    std::vector<double> currPoint;
    const data_t& inputData = input.getData();

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

Points2D& JarvisScan2D::solvePara(const Points2D& input, Points2D& output)
{
    if (input.getSize() <= 1) {
        output = input;
        return output;
    }

    const data_t& inputData = input.getData();

    // find extremes
    unsigned maxXIdx = 0, maxYIdx = 0, minXIdx = 0, minYIdx = 0;
    for (unsigned i = 1; i < inputData.size(); i++) {
        double minXdif = inputData[i][0] - inputData[minXIdx][0],
               minYdif = inputData[i][1] - inputData[minYIdx][1],
               maxXdif = inputData[i][0] - inputData[maxXIdx][0],
               maxYdif = inputData[i][1] - inputData[maxYIdx][1];
        if (fabs(minXdif) < EPS) {
            if (inputData[i][1] < inputData[minXIdx][1]) {
                minXIdx = i;
            }
        } else if (minXdif < -EPS) {
            minXIdx = i;
        }
        if (fabs(minYdif) < EPS) {
            if (inputData[i][0] < inputData[minYIdx][0]) {
                minYIdx = i;
            }
        } else if (minYdif < -EPS) {
            minYIdx = i;
        }
        if (fabs(maxXdif) < EPS) {
            if (inputData[i][1] > inputData[maxXIdx][1]) {
                maxXIdx = i;
            }
        } else if (maxXdif > EPS) {
            maxXIdx = i;
        }
        if (fabs(maxYdif) < EPS) {
            if (inputData[i][0] > inputData[maxYIdx][0]) {
                maxYIdx = i;
            }
        } else if (maxYdif > EPS) {
            maxYIdx = i;
        }
    }

    data_t part[4];

#pragma omp parallel
    {
#pragma omp sections
        {
#pragma omp section
            scan(inputData, part[0], minXIdx, maxYIdx);
#pragma omp section
            scan(inputData, part[1], maxYIdx, maxXIdx);
#pragma omp section
            scan(inputData, part[2], maxXIdx, minYIdx);
#pragma omp section
            scan(inputData, part[3], minYIdx, minXIdx);
        }
    }

    for (int i = 0; i < 4; i++) {
        for (auto& pt : part[i]) {
            output.add(pt);
        }
    }

    if (output.getSize() == 0) {
        output.add(inputData[0]);
    }

    return output;
}

void JarvisScan2D::scan(const data_t& input, data_t& output,
                        unsigned beginIdx, unsigned endIdx)
{
    unsigned currIdx = beginIdx, nextIdx;
    while (currIdx != endIdx) {
        output.push_back(input[currIdx]);
        // avoid setting same point as next
        nextIdx = !currIdx;

        // check orientation for all remaining n - 1 points
        for (unsigned i = 0; i < input.size(); i++) {
            if (i == currIdx) {
                continue;
            }

            int o = orientation(input[currIdx][0],
                                input[currIdx][1],
                                input[nextIdx][0],
                                input[nextIdx][1],
                                input[i      ][0],
                                input[i      ][1]);

            if (o == 0) {
                // exclude collinear points
                if (dist(input[currIdx][0],
                         input[currIdx][1],
                         input[i][0],
                         input[i][1])
                    >
                    dist(input[currIdx][0],
                         input[currIdx][1],
                         input[nextIdx][0],
                         input[nextIdx][1])) {
                    nextIdx = i;
                }
            } else if (o == 1) {
                // point to the left of current hull face
                nextIdx = i;
            }
        }
        currIdx = nextIdx;
    }
}

}
