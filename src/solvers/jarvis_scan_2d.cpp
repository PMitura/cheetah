#include "jarvis_scan_2d.h"

namespace ch {

Points2D& JarvisScan2D::solve(const Points2D& input, Points2D& output)
{
    std::vector<double> nextPoint;
    const data_t& inputData = input.getData();

    // find min X
    int minIndex = 0;
    for (int i = 1; i < inputData.size(); i++) {
        if (inputData[i][0] < inputData[minIndex][0]) {
            minIndex = i;
        }
    }
    nextPoint = inputData[minIndex];

    return output;
}

}
