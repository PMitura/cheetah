#include "approximators/bfp2d.h"

namespace ch
{

BFP2D::BFP2D()
{
    name_ = "BFP";
}

Points2D& BFP2D::approximate(const Points2D& input, Points2D& output)
{
    return sequential(input, output);
}


Points2D& BFP2D::sequential(const Points2D& input, Points2D& output)
{
    if (input.getSize() <= 2) {
        output = input;
        return output;
    }

    // fixed parameter, higher number means more accuracy
    const int stripsCount = 100;
    const data_t& inputData = input.getData();

    // find extremes
    unsigned minXminY = 0, minXmaxY = 0,
             maxXminY = 0, maxXmaxY = 0;

    double xD;
    for (unsigned i = 1; i < inputData.size(); i++) {
        xD = inputData[i][0] - inputData[minXminY][0];
        if (xD < EPS) {
            if (xD < -EPS) {
                minXminY = i;
                minXmaxY = i;
            } else {
                if (inputData[i][1] < inputData[minXminY][1] - EPS) {
                    minXminY = i;
                } else if (inputData[i][1] > inputData[minXmaxY][1] + EPS) {
                    minXmaxY = i;
                }
            }
        }

        xD = inputData[i][0] - inputData[maxXminY][0];
        if (xD > -EPS) {
            if (xD > EPS) {
                maxXminY = i;
                maxXmaxY = i;
            } else {
                if (inputData[i][1] < inputData[maxXminY][1] - EPS) {
                    maxXminY = i;
                } else if (inputData[i][1] > inputData[maxXmaxY][1] + EPS) {
                    maxXmaxY = i;
                }
            }
        }
    }
    double minX = inputData[minXminY][0],
           maxX = inputData[maxXminY][0];

    // minX = maxX
    if (fabs(minX - maxX) < EPS) {
        output.add(inputData[minXmaxY]);
        if (fabs(inputData[minXmaxY][1] - inputData[minXminY][1]) > EPS) {
            output.add(inputData[minXminY]);
        }
        return output;
    }

    // indexes of min/max in strips
    std::vector<std::pair<int, int>> strips(stripsCount + 2);
    strips.front() = {minXminY, minXmaxY};
    strips.back()  = {maxXminY, maxXmaxY};
    for (int i = 1; i < stripsCount - 1; i++) {
        strips[i] = {-1, -1};
    }

    // sort points into strips
    unsigned sIdx;
    for (int i = 0; i < inputData.size(); i++) {
        const point_t& curr = inputData[i];
        // first / last strip
        if (   fabs(curr[0] - inputData[minXminY][0]) < EPS
            || fabs(curr[0] - inputData[maxXminY][0]) < EPS) {
            continue;
        }
        if (!ccw(inputData[minXminY], inputData[maxXminY], curr)) {
            // below low
            sIdx = stripsCount * (curr[0] - minX) / (maxX - minX) + 1;
            if (strips[sIdx].first == -1) {
                strips[sIdx].first = i;
            } else if (curr[1] > inputData[strips[sIdx].first][1] + EPS) {
                strips[sIdx].first = i;
            }
        } else if (ccw(inputData[minXmaxY], inputData[maxXmaxY], curr)) {
            // above high
            sIdx = stripsCount * (curr[0] - minX) / (maxX - minX) + 1;
            if (strips[sIdx].second == -1) {
                strips[sIdx].second = i;
            } else if (curr[1] > inputData[strips[sIdx].second][1] + EPS) {
                strips[sIdx].second = i;
            }
        }
    }

    // scan lower - same as in monochain
    std::vector<unsigned> pStack(stripsCount + 3);
    int sSize = 0, curr;
    for (int i = 0; i < stripsCount + 2; i++) {
        curr = strips[i].first;
        if (curr == -1) {
            continue;
        }
        while (sSize >= 2 && !ccw(inputData[pStack[sSize - 2]],
                                  inputData[pStack[sSize - 1]],
                                  inputData[curr])) {
            sSize--;
        }
        pStack[sSize++] = curr;
    }
    for (int i = 0; i < sSize; i++) {
        output.add(inputData[pStack[i]]);
    }

    // scan upper
    sSize = 0; pStack.clear();
    for (int i = stripsCount + 1; i >= 0; i++) {
        curr = strips[i].second;
        if (curr == -1) {
            continue;
        }
        while (sSize >= 2 && !ccw(inputData[pStack[sSize - 2]],
                                  inputData[pStack[sSize - 1]],
                                  inputData[curr])) {
            sSize--;
        }
        pStack[sSize++] = curr;
    }
    for (int i = 0; i < sSize; i++) {
        output.add(inputData[pStack[i]]);
    }

    return output;
}

}
