#include "solvers/graham_scan_2d.h"

namespace ch
{

GrahamScan2D::GrahamScan2D()
{
    name_ = "Graham Scan";
}

Points2D& GrahamScan2D::solve(const Points2D& input, Points2D& output)
{
    // return solveSequential(input, output);
    return solveParallel(input, output);
}

Points2D& GrahamScan2D::solveSequential(const Points2D& input,
                                        Points2D& output)
{
    if (input.getSize() <= 2) {
        output = input;
        return output;
    }

    order_.clear(); polar_.clear(); output.clear();
    const data_t& inputData = input.getData();

    for (unsigned i = 0; i < inputData.size(); i++)
        order_.push_back(i);
    pivot_ = findMinY(inputData);
    std::swap(order_.at(0), order_.at(pivot_));

    sortPoints(inputData);

    scan(inputData, output);

    return output;
}

Points2D& GrahamScan2D::solveParallel(const Points2D& input, Points2D& output)
{
    if (input.getSize() <= 2) {
        output = input;
        return output;
    }

    order_.clear(); polar_.clear(); output.clear();
    const data_t& inputData = input.getData();

    for (unsigned i = 0; i < inputData.size(); i++)
        order_.push_back(i);
    pivot_ = findMinY(inputData);
    std::swap(order_[0], order_[pivot_]);

    sortPointsParallel(inputData);

    scan(inputData, output);

    return output;

}

int GrahamScan2D::findMinY(const data_t& points)
{
    int minIndex = 0;
    for (unsigned i = 1; i < points.size(); i++) {
        double delta = points[minIndex][1] - points[i][1];
        if (delta > EPS) {
            minIndex = i;
        } else if (fabs(delta) < EPS) {
            if (points[minIndex][0] < points[i][0]) {
                minIndex = i;
            }
        }
    }
    return minIndex;
}

void GrahamScan2D::computeAngles(const data_t& points)
{
    double dx, dy;
    for (unsigned i = 0; i < points.size(); i++) {
        dx = points[i][0] - points[pivot_][0];
        dy = points[i][1] - points[pivot_][1];
        polar_.push_back(-(dx / dy));
    }
}

void GrahamScan2D::scan(const data_t& inputData, Points2D& output)
{
    unsigned iPtr = 2, sPtr = 2, iSize = inputData.size();
    unsigned * ptStack = new unsigned[iSize];
    ptStack[0] = order_[0];
    ptStack[1] = order_[1];
    while (iPtr < iSize) {
        if (ccw(inputData[ptStack[0]][0],
                inputData[ptStack[0]][1],
                inputData[ptStack[1]][0],
                inputData[ptStack[1]][1],
                inputData[order_[iPtr]][0],
                inputData[order_[iPtr]][1])) {
            break;
        }
        ptStack[1] = order_[iPtr++];
    }

    while (iPtr < iSize) {
        if (ccw(inputData[ptStack[sPtr - 2]][0],
                inputData[ptStack[sPtr - 2]][1],
                inputData[ptStack[sPtr - 1]][0],
                inputData[ptStack[sPtr - 1]][1],
                inputData[order_[iPtr]][0],
                inputData[order_[iPtr]][1])) {
            ptStack[sPtr++] = order_[iPtr++];
        } else {
            sPtr--;
        }
    }

    // handle last point collinear with first
    if (sPtr >= 3) {
        if (ccw(inputData[ptStack[sPtr - 2]][0],
                inputData[ptStack[sPtr - 2]][1],
                inputData[ptStack[sPtr - 1]][0],
                inputData[ptStack[sPtr - 1]][1],
                inputData[order_[0]][0],
                inputData[order_[0]][1]) == 0) {
            sPtr--;
        }
    }

    for (unsigned i = 0; i < sPtr; i++) {
        output.add(inputData[ptStack[i]]);
    }
    delete[] ptStack;
}

void GrahamScan2D::sortPoints(const data_t& inputData)
{
    computeAngles(inputData);
    std::stable_sort((order_.begin()) + 1, order_.end(),
            AngleCmp(*this, inputData));
}

void GrahamScan2D::sortPointsParallel(const data_t& inputData)
{
    computeAngles(inputData);
    __gnu_parallel::stable_sort((order_.begin()) + 1, order_.end(),
            AngleCmp(*this, inputData));
}

bool GrahamScan2D::AngleCmp::operator()(const unsigned& a, const unsigned& b)
{
    double x = part_.polar_[a] - part_.polar_[b];
    if (fabs(x) < EPS) {
        return   dist({data_[part_.pivot_][0], data_[part_.pivot_][1]},
                        {data_[a][0], data_[a][1]})
               < dist({data_[part_.pivot_][0], data_[part_.pivot_][1]},
                        {data_[b][0], data_[b][1]});
    }
    return x < EPS;
}


}
