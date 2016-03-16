#include "solvers/graham_scan_2d.h"

namespace ch
{

GrahamScan2D::GrahamScan2D()
{
    order_ = new std::vector<int>;
    polar_ = new std::vector<double>;
}

GrahamScan2D::~GrahamScan2D()
{
    delete order_;
    delete polar_;
}

Points2D& GrahamScan2D::solve(const Points2D& input, Points2D& output)
{
    order_ -> clear();
    polar_ -> clear();
    const data_t& inputData = input.getData();

    for (int i = 0; i < input.getSize(); i++)
        order_ -> push_back(i);
    int pivot = findMinY(inputData);
    std::swap(order_[0], order_[pivot]);

    computeAngles(pivot);
}

int GrahamScan2D::findMinY(const data_t& points)
{
    int minIndex = 0;
    for (int i = 1; i < points.size(); i++) {
        if (points[i][1] < points[minIndex][1]) {
            minIndex = i;
        }
    }
    return minIndex;
}

void GrahamScan2D::computeAngles(int pivotIdx)
{

}

}
