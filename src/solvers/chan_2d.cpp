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
    // find hull size
    for (unsigned h = 1; ppow(h) < input.getSize(); h++) {
        std::vector<Points2D> hulls;
        findHulls(input, hulls, ppow(h));

        unsigned pivot,
                 minHull = findMinHull(hulls, pivot),
                 candidate = (pivot + 1) % hulls[minHull].getSize();
        point_t ppt = hulls[minHull].getData()[pivot],
                cpt = hulls[minHull].getData()[candidate];
        std::vector<std::pair<unsigned, unsigned>> overallHull;
        // iterate at most h times (guessed number of hull points)
        for (unsigned i = 0; i < hulls.size(); i++) {
            findNext
        }
    }
    return output;
}

void Chan2D::findHulls(const Points2D& input, std::vector<Points2D>& hulls,
                       unsigned step)
{
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

std::pair<unsigned, unsigned> Chan2D::findNext(std::vector<Points2D>& hulls,
        std::pair<unsigned, unsigned> curr)
{
    
}

unsigned Chan2D::findMinHull(std::vector<Points2D>& hulls, unsigned& minPt)
{
    unsigned lies = 0;
    double minX = DBL_MIN;
    for (unsigned h = 0; h < hulls.size(); h++) {
        const data_t& hull = hulls[h].getData();
        for (unsigned i = 0; i < hull.size(); i++) {
            if (hull[i][0] < minX) {
                minX = hull[i][0];
                minPt = i;
                lies = h;
            }
        }
    }
    return lies;
}

unsigned Chan2D::findTangent(const Points2D& hull, point_t p)
{
    return 0;
}

}
