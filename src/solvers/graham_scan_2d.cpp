#include "solvers/graham_scan_2d.h"

namespace ch
{

GrahamScan2D::GrahamScan2D()
{
    order_ = new std::vector<unsigned>;
    polar_ = new std::vector<double>;
}

GrahamScan2D::~GrahamScan2D()
{
    delete order_;
    delete polar_;
}

Points2D& GrahamScan2D::solve(const Points2D& input, Points2D& output)
{
    if (input.getSize() <= 2) {
        output = input;
        return output;
    }

    order_ -> clear();
    polar_ -> clear();
    output.clear();
    const data_t& inputData = input.getData();

    for (unsigned i = 0; i < inputData.size(); i++)
        order_ -> push_back(i);
    pivot_ = findMinY(inputData);
    std::swap(order_ -> at(0), order_ -> at(pivot_));
    computeAngles(inputData);

    std::sort((order_ -> begin()) + 1, order_ -> end(), 
            AngleCmp(*this, inputData));

    std::vector<unsigned>* ptStack = new std::vector<unsigned>;
    ptStack -> push_back(order_ -> at(0));
    ptStack -> push_back(order_ -> at(1));
    unsigned iPtr = 2, sPtr = 2, iIdx, sIdx1, sIdx2;
    while (iPtr < inputData.size()) {
        iIdx  = order_  -> at(iPtr);
        sIdx1 = ptStack -> at(1);
        sIdx2 = ptStack -> at(0);
        if (orientation(inputData[sIdx2][0], inputData[sIdx2][1],
                        inputData[sIdx1][0], inputData[sIdx1][1],
                        inputData[iIdx ][0], inputData[iIdx ][1]) == 1) {
            break;
        }
        iPtr++;
        ptStack -> at(1) = iIdx;
    }
    while (iPtr < inputData.size()) {
        iIdx  = order_  -> at(iPtr);
        sIdx1 = ptStack -> at(sPtr - 1);
        sIdx2 = ptStack -> at(sPtr - 2);
        if (orientation(inputData[sIdx2][0], inputData[sIdx2][1],
                        inputData[sIdx1][0], inputData[sIdx1][1],
                        inputData[iIdx ][0], inputData[iIdx ][1]) == 1) {
            ptStack -> push_back(iIdx);
            sPtr++;
            iPtr++;
        } else {
            ptStack -> pop_back();
            sPtr--;
        }
    }

    // handle last point collinear with first
    if (ptStack -> size() >= 3) {
        unsigned iIdx = order_  -> at(0),
                sIdx1 = ptStack -> at(sPtr - 1),
                sIdx2 = ptStack -> at(sPtr - 2);
        if (orientation(inputData[sIdx2][0], inputData[sIdx2][1],
                        inputData[sIdx1][0], inputData[sIdx1][1],
                        inputData[iIdx ][0], inputData[iIdx ][1]) != 1) {
            ptStack -> pop_back();
        }
    }

    for (auto& i : (*ptStack)) {
        output.add(inputData[i]);
    }

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
        polar_ -> push_back(-(dx / dy));
    }
}

bool GrahamScan2D::AngleCmp::operator()(const int& a, const int& b)
{
    double x = (*part_.polar_)[a] - (*part_.polar_)[b];
    if (fabs(x) < EPS) {
        return   dist({data_[part_.pivot_][0], data_[part_.pivot_][1]},
                        {data_[a][0], data_[a][1]})
               < dist({data_[part_.pivot_][0], data_[part_.pivot_][1]},
                        {data_[b][0], data_[b][1]});
    }
    return x < EPS;
}


}
