#include "solvers/jarvis_scan_3d.h"

namespace ch
{


Polyhedron& JarvisScan3D::solve(const Points3D& input, Polyhedron& output)
{
    // edge cases
    if (input.getSize() == 0) {
        return output;
    } else if (input.getSize() < 3) {
        Points2D plane;
        const data_t& inputData = input.getData();
        for (auto& i : inputData) {
            plane.add(i);
        }
        output.addFace(plane);
        return output;
    }

    solveNaive(input, output);

    return output;
}

Polyhedron& JarvisScan3D::solveNaive(const Points3D& input, Polyhedron& output)
{
    const data_t& idata = input.getData();
    std::pair<unsigned, unsigned> init = findInitial(idata);
    if (init.first >= init.second) {
        std::swap(init.first, init.second);
    }
    // unprocessed edges
    std::set<std::pair<unsigned, unsigned>> edges;
    edges.insert(init);

    while (!edges.empty()) {
        // points a, b
        std::pair<unsigned, unsigned> curr = *(edges.begin());
        point_t vab = {idata[curr.first][0] - idata[curr.second][0],
                       idata[curr.first][1] - idata[curr.second][1],
                       idata[curr.first][2] - idata[curr.second][2]};

        point_t vcb, vperp;
        unsigned c = -1;
        // find some non-collinear c
        for (unsigned i = 0; i < idata.size(); i++) {
            // cannot reuse a, b
            if (i == curr.first || i == curr.second) {
                continue;
            }
            vcb = {idata[i][0] - idata[curr.second][0],
                   idata[i][1] - idata[curr.second][1],
                   idata[i][2] - idata[curr.second][2]};
            // find vector perpendicular to abc
            vperp = perpend3d(vab, vcb);
            // if c not collinear, use it
            if (!zeroVect(vperp)) {
                c = i;
                break;
            }
        }

        // all input points on one line
        if (c == -1) {
            return output;
        }

        unsigned prevC, currC = c, itN = 0;
        std::vector<unsigned> onPlane;
        // find such c, that all points are on one side of plane abc
        do {
            // update perpendicular vector
            vperp = perpendNormal3d(vab, vcb);
            // find new c
            double dDist, maxDist;
            onPlane.clear();
            prevC = currC;
            for (unsigned i = 0; i < idata.size(); i++) {
                // cannot reuse a, b, c
                if (i == curr.first || i == curr.second || i == c) {
                    continue;
                }

                // compute directional distance
                dDist = dot(idata[i][0] - idata[curr.second][0],
                            idata[i][1] - idata[curr.second][1],
                            idata[i][2] - idata[curr.second][2],
                            vperp[0], vperp[1], vperp[2]);
                if (dDist > maxDist) {
                    currC = i;
                    maxDist = dDist;
                } else if (fabs(dDist) < EPS) {
                    // on the same plane as abc
                    onPlane.push_back(i);
                }
            }

            // detect infinite loop
            if (itN++ >= idata.size()) {
                // this happening means wrong choice of ab
                // report occurencies of this call as a bug
                std::cerr << "Cannot find third point of plane" << std::endl;
                return output;
            }

            vcb = {idata[c][0] - idata[curr.second][0],
                   idata[c][1] - idata[curr.second][1],
                   idata[c][2] - idata[curr.second][2]};
            // find new perpendicular vector
        } while (prevC != currC); // c not changed => all points on one side

        // PROBLEM
        // this returns all points on faces, but I only need their convex hull
        // we need a way to find convex hull od points on plane, but in 3d
    }

    return output;
}

double JarvisScan3D::randomOne()
{
    return (rand() % 2) ? -1 : 1;
}

std::pair<unsigned, unsigned> JarvisScan3D::findInitial(const data_t& input)
{
    srand(time(NULL));
    int rndx, rndy, rndz;
    double maxd;
    unsigned far, farcnt;

    // find farthest point in some direction
    do {
        far = 0;
        farcnt = 1;
        maxd = DBL_MIN;
        // random direction vector
        rndx = randomOne(); rndy = randomOne(); rndz = randomOne();
        // R(rndx << ", " << rndy << ", " << rndz);
        for (unsigned i = 0; i < input.size(); i++) {
            double d = dot(input[i][0], input[i][1], input[i][2], 
                           rndx,        rndy,        rndz);
            double dif = d - maxd;
            if (dif > EPS) {
                maxd = d;
                farcnt = 1;
                far = i;
            } else if (fabs(dif) < EPS) {
                farcnt++;
            }
        }
    } while (farcnt > 3);

    unsigned paired = 0;
    maxd = DBL_MAX;
    double minLen = DBL_MAX;
    // find second point on edge = point with smallest angle to plane
    // perpendicular to rng vector from previous step
    for (unsigned i = 0; i < input.size(); i++) {
        if (i == far) {
            continue;
        }
        point_t vFarI = {input[far][0] - input[i][0],
                         input[far][1] - input[i][1],
                         input[far][2] - input[i][2]};
        double iLen = sqrt(vFarI[0]*vFarI[0] + vFarI[1]*vFarI[1]
                           + vFarI[2]*vFarI[2]);
        // directional distance
        double d = dot(vFarI[0], vFarI[1], vFarI[2], 
                       rndx,     rndy,     rndz);
        d /= iLen;
        double dif = d - maxd;
        if (dif < -EPS) {
            paired = i;
            maxd = d;
            minLen = iLen;
        } else if (fabs(dif) < EPS) {
            if (iLen < minLen - EPS) {
                paired = i;
                minLen = iLen;
            }
        }
    }

    return {far, paired};
}

}
