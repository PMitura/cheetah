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
    if (input.getSize() <= 2) {
        output = input;
        return output;
    }

    // find hull size, stop one step after input size
    for (unsigned h = 1; ppow(h-1) < input.getSize(); h++) {
        std::vector<Points2D> hulls;
        findHulls(input, hulls, ppow(h));

        unsigned pivot = 0,
                 minHull = findMinHull(hulls, pivot);
        std::vector<std::pair<unsigned, unsigned>> overallHull;
        std::pair<unsigned, unsigned> curr = {pivot, minHull};
        overallHull.push_back(curr);

        // iterate at most 2^2^h times (guessed number of hull points)
        bool fnd = 0;
        for (unsigned i = 0; i < ppow(h); i++) {
            curr = findNext(hulls, curr);
            if (curr == overallHull[0]) {
                fnd = 1;
                break;
            }
            overallHull.push_back(curr);
        }
        if (fnd) {
            for (auto& i : overallHull) {
                output.add(hulls[i.second].getData()[i.first]);
            }
            return output;
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
            part.add(inputData[j]);
        }
        hulls.push_back(Points2D());
        solver.solve(part, hulls.back());
    }
}

std::pair<unsigned, unsigned> Chan2D::findNext(std::vector<Points2D>& hulls,
        std::pair<unsigned, unsigned> curr)
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
            double dif = hull[i][0] - minX;
            if (dif < -EPS) {
                minX = hull[i][0];
                minPt = i;
                lies = h;
            } else if (fabs(dif) < EPS) {
                if (hull[i][1] < hulls[lies].getData()[minPt][1]) {
                    minX = hull[i][0];
                    minPt = i;
                    lies = h;
                }
            }
        }
    }
    return lies;
}

unsigned Chan2D::findTangent(const Points2D& hull, point_t& p)
{
    const data_t& hdata = hull.getData();
    if (hdata.size() == 1) {
        return 0;
    }

    if (hdata.size() == 2) {
        if (dist(p, hdata[0]) > dist(p, hdata[1])) {
            return 0;
        } else {
            return 1;
        }
    }

    unsigned left = 0, right = hdata.size(), mid, s = hdata.size();
    int olFrnt = orientation(p[0],            p[1],
                             hdata[0][0],     hdata[0][1],
                             hdata[1][0],     hdata[1][1]),
        olBack = orientation(p[0],            p[1],
                             hdata[0][0],     hdata[0][1],
                             hdata[s - 1][0], hdata[s - 1][1]);

    while (left < right) {
        mid = (left + right) / 2;
        // R(left << " - " << mid << " - " << right);
        int dm = (mid == 0) ? s - 1 : mid - 1;

        int omFrnt = orientation(p[0],                  p[1],
                                 hdata[mid][0],         hdata[mid][1],
                                 hdata[(mid+1) % s][0], hdata[(mid+1) % s][1]);
        int omBack = orientation(p[0],          p[1],
                                 hdata[mid][0], hdata[mid][1],
                                 hdata[dm][0],  hdata[dm][1]);
        int omSelf = orientation(p[0],           p[1],
                                 hdata[left][0], hdata[left][1],
                                 hdata[mid][0],  hdata[mid][1]);

        if (omBack == 1 && omFrnt == 1) {
            return mid;
        }

        if (   (omSelf == 2 && omBack == 2)
            || (omSelf == 1 && (olFrnt == 2 || olBack == olFrnt))) {
            right = mid;
        } else {
            left = mid + 1;
            olBack = omFrnt;
            if (omFrnt == 1) {
                olBack = 2;
            } else if (omFrnt == 2) {
                olBack = 1;
            }
            olFrnt = orientation(p[0], p[1],
                hdata[left % s][0], hdata[left % s][1],
                hdata[(left+1) % s][0], hdata[(left+1) % s][1]);
        }
    }
    return 1;
}

}
