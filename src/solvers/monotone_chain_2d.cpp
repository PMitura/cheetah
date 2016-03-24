#include "solvers/monotone_chain_2d.h"

namespace ch
{

MonotoneChain2D::MonotoneChain2D()
{
    name_ = "Monotone Chain";
}

Points2D& MonotoneChain2D::solve(const Points2D& input, Points2D& output)
{
    return solveSequential(input, output);
}

Points2D& MonotoneChain2D::solveSequential(const Points2D& input,
                                           Points2D& output)
{
    if (input.getSize() <= 1) {
        output = input;
        return output;
    }

    data_t inputData = input.getData();
    double tA = omp_get_wtime();
    std::sort(inputData.begin(), inputData.end(), PointCmp());
    double tB = omp_get_wtime();
    R("") D(tB - tA);

    unsigned * lower = new unsigned[input.getSize()],
             * upper = new unsigned[input.getSize()];

    unsigned lowerSize = scanLower(inputData, lower);
    unsigned upperSize = scanUpper(inputData, upper);

    // do not include last points to avoid duplicates
    for (unsigned i = 0; i < lowerSize - 1; i++) {
        output.add(inputData[lower[i]]);
    }
    delete[] lower;

    for (unsigned i = 0; i < upperSize - 1; i++) {
        output.add(inputData[upper[i]]);
    }
    delete[] upper;

    return output;
}

Points2D& MonotoneChain2D::solveParallel(const Points2D& input,
                                         Points2D& output)
{

    return output;
}

unsigned MonotoneChain2D::scanLower(const data_t& input, unsigned* lower)
{
    unsigned sSize = 0;

    for (unsigned i = 0; i < input.size(); i++) {
        while (sSize >= 2 && !ccw(input[lower[sSize - 2]][0],
                                  input[lower[sSize - 2]][1],
                                  input[lower[sSize - 1]][0],
                                  input[lower[sSize - 1]][1],
                                  input[i][0],
                                  input[i][1])) {
            sSize--;
        }
        lower[sSize++] = i;
    }

    return sSize;
}

unsigned MonotoneChain2D::scanUpper(const data_t& input, unsigned* upper)
{
    unsigned sSize = 0;

    for (int i = (int) input.size() - 1; i >= 0; i--) {
        while (sSize >= 2 && !ccw(input[upper[sSize - 2]][0],
                                  input[upper[sSize - 2]][1],
                                  input[upper[sSize - 1]][0],
                                  input[upper[sSize - 1]][1],
                                  input[i][0],
                                  input[i][1])) {
            sSize--;
        }
        upper[sSize++] = i;
    }

    return sSize;
}

bool MonotoneChain2D::PointCmp::operator()(const point_t& a,
                                           const point_t& b)
{
    double dif = a[0] - b[0];
    if (fabs(dif) < EPS) {
        return a[1] < b[1];
    }
    return dif > EPS;
}

}
