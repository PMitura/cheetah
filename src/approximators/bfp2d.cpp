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

    const data_t& inputData = input.getData();

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

    // minX = maxX
    if (fabs(inputData[minXmaxY][0] - inputData[maxXmaxY][0]) < EPS) {
        output.add(inputData[minXmaxY]);
        if (fabs(inputData[minXmaxY][1] - inputData[minXminY][1]) > EPS) {
            output.add(inputData[minXminY]);
        }
        return output;
    }



    return output;
}

}
