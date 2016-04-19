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

        unsigned pivot = 0,
                 minHull = findMinHull(hulls, pivot);
        std::vector<std::pair<unsigned, unsigned>> overallHull;
        std::pair<unsigned, unsigned> curr = {pivot, minHull};
        overallHull.push_back(curr);

        // iterate at most h times (guessed number of hull points)
        for (unsigned i = 0; i < hulls.size(); i++) {
            curr = findNext(hulls, curr);

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
        std::pair<unsigned, unsigned>& curr)
{
    std::pair<unsigned, unsigned> cand;
    cand.first = (curr.first + 1) % (hulls[curr.second].getSize());
    cand.second = curr.second;
    point_t currP = hulls[curr.second].getData()[curr.first],
            candP = hulls[cand.second].getData()[cand.first],
            propP;

    for (unsigned sub = 0; sub < hulls.size(); sub++) {
        if (sub == curr.second) {
            continue;
        }

        unsigned tgt = findTangent(hulls[sub], currP);
        propP = hulls[sub].getData()[tgt];
        int o = orientation(currP[0], currP[1],
                            candP[0], candP[1],
                            propP[0], propP[1]);

        if (o == 2) {
            // right turn
            cand = {tgt, sub};
            candP = propP;
        } else if (o == 0) {
            // collinear
            if (dist(currP, candP) < dist(currP, propP)) {
                cand = {tgt, sub};
                candP = propP;
            }
        }
    }

    return cand;
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

unsigned Chan2D::findTangent(const Points2D& hull, point_t& p)
{
    return 0;
}

}
