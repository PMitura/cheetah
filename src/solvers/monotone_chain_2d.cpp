#include "solvers/monotone_chain_2d.h"

namespace ch
{

MonotoneChain2D::MonotoneChain2D()
{
    name_ = "Monotone Chain";
}

Points2D& MonotoneChain2D::solve(const Points2D& input, Points2D& output)
{
    // return solveSequential(input, output);
    return solveParallel(input, output);
}

Points2D& MonotoneChain2D::solveSequential(const Points2D& input,
                                           Points2D& output)
{
    if (input.getSize() <= 2) {
        output = input;
        return output;
    }

    const data_t& inputData = input.getData();
    order_.clear();
    for (unsigned i = 0; i < inputData.size(); i++) {
        order_.push_back(i);
    }

    sortPtsDirect(inputData);

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
    if (input.getSize() <= 2) {
        output = input;
        return output;
    }

    const data_t& inputData = input.getData();
    order_.clear();
    for (unsigned i = 0; i < inputData.size(); i++) {
        order_.push_back(i);
    }

    // double tA = omp_get_wtime();
    sortPtsParallel(inputData);
    // double tB = omp_get_wtime();
    // R("") R("sort time:  " << tB - tA << " ms") std::cout << "  total time: ";

    unsigned * lower = new unsigned[input.getSize()],
             * upper = new unsigned[input.getSize()];

    unsigned lowerSize, upperSize;
#pragma omp parallel
    {
#pragma omp sections
        {
#pragma omp section
            lowerSize = scanLower(inputData, lower);
#pragma omp section
            upperSize = scanUpper(inputData, upper);
        }
    }

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
                                  input[order_[i]][0],
                                  input[order_[i]][1])) {
            sSize--;
        }
        lower[sSize++] = order_[i];
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
                                  input[order_[i]][0],
                                  input[order_[i]][1])) {
            sSize--;
        }
        upper[sSize++] = order_[i];
    }

    return sSize;
}

void MonotoneChain2D::sortPtsDirect(const data_t& input)
{
    std::sort(order_.begin(), order_.end(), PointCmpDirect(input));
}

void MonotoneChain2D::sortPtsCache(const data_t& input)
{
    xar_.clear(); yar_.clear();
    for (unsigned i = 0; i < input.size(); i++) {
        xar_.push_back(input[i][0]);
        yar_.push_back(input[i][1]);
    }
    std::sort(order_.begin(), order_.end(), PointCmpCache(*this, input));
}

void MonotoneChain2D::sortPtsParallel(const data_t& input)
{
    xar_.clear(); yar_.clear();
    for (unsigned i = 0; i < input.size(); i++) {
        xar_.push_back(input[i][0]);
        yar_.push_back(input[i][1]);
    }
    __gnu_parallel::stable_sort(order_.begin(), order_.end(),
                                PointCmpCache(*this, input));
}

bool MonotoneChain2D::PointCmpDirect::operator()(const unsigned& a,
                                                const unsigned& b)
{
    double dif = data_[a][0] - data_[b][0];
    if (fabs(dif) < EPS) {
        return data_[a][1] > data_[b][1];
    }
    return dif > EPS;
}

bool MonotoneChain2D::PointCmpCache::operator()(const unsigned& a,
                                                const unsigned& b)
{
    double dif = part_.xar_[a] - part_.xar_[b];
    if (fabs(dif) < EPS) {
        return part_.yar_[a] > part_.yar_[b];
    }
    return dif > EPS;
}

}
