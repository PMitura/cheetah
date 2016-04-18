#include "solvers/chan_2d.h"

namespace ch
{

Chan2D::Chan2D()
{
    name_ = "Chan";
}

Points2D& Chan2D::solve(const Points2D& input, Points2D& output)
{
    return solveNaive(input, output);
}

Points2D& Chan2D::solveNaive(const Points2D& input, Points2D& output)
{
    for (unsigned h = 1; (1U << (1U << h)) < input.getSize(); h++) {
        std::vector<Points2D> hulls;
        findHulls(input, hulls, (1U << (1U << h)));
    }
    return output;
}

void Chan2D::findHulls(const Points2D& input, std::vector<Points2D>& hulls,
                       unsigned hullSize)
{
    unsigned step = input.getSize() / hullSize;
    const data_t& inputData = input.getData();
    for (unsigned i = 0; i < input.getSize(); i += step) {
        Points2D part;
        GrahamScan2D solver;
        for (unsigned j = i; j < std::min(i + step, input.getSize()); j++) {
            part.add(inputData[i]);
        }
        hulls.push_back(Points2D());
        solver.solve(part, hulls.back());
    }
}

unsigned findMinHull(std::vector<Points2D>& hulls, point_t& minPt)
{
    unsigned lies = 0;
    double minX = DBL_MIN;
    for (unsigned h = 0; h < hulls.size(); h++) {
        const data_t& hull = hulls[h].getData();
        for (auto& pt : hull) {
            if (pt[0] < minX) {
                minX = pt[0];
                minPt = pt;
                lies = h;
            }
        }
    }
    return lies;
}

}
