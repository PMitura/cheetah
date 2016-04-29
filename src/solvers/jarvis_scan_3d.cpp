#include "solvers/jarvis_scan_3d.h"

namespace ch
{

JarvisScan3D::JarvisScan3D()
{
    // fixed value is the lesser evil here, relative epsilon causes all sorts
    // of weird behavior
    EPS_LOC = 1e-6;
}

Polyhedron& JarvisScan3D::solve(const Points3D& input, Polyhedron& output)
{
    // edge cases
    if (input.getSize() == 0) {
        return output;
    } else if (input.getSize() < 3) {
        Points3D plane;
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
    typedef std::pair<unsigned, unsigned> edge_t;
    edge_t init = findInitial(idata);
    // discovered and processed edges
    std::set<edge_t> fresh, closed;

    do {
        // points a, b
        edge_t curr = {init.first, init.second};
        if (fresh.size() > 0) {
            curr = *(fresh.begin());
        }
        point_t vab = {idata[curr.first][0] - idata[curr.second][0],
                       idata[curr.first][1] - idata[curr.second][1],
                       idata[curr.first][2] - idata[curr.second][2]};

        R(" ");
        R("OVER EDGE: " << idata[curr.first][0] << ", " <<
                           idata[curr.first][1] << ", " <<
                           idata[curr.first][2] << " | " <<
                           idata[curr.second][0] << ", " <<
                           idata[curr.second][1] << ", " <<
                           idata[curr.second][2] << 
                           " (id: " << curr.first << ", " << curr.second << ")");

        point_t vcb, vperp;
        unsigned c = UINT_MAX;
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
        if (c == UINT_MAX) {
            R("collinear set");
            return output;
        }


        unsigned prevC, currC = c, itN = 0;
        std::vector<unsigned> onPlane;
        // find such c, that all points are on one side of plane abc
        do {
            // update perpendicular vector
            vperp = perpendNormal3d(vab, vcb);
            // find new c
            double dDist, maxDist = 0;
            onPlane.clear();
            onPlane.push_back(currC);
            prevC = currC;
            R("TRY C: " << idata[currC][0] << ", " <<
                           idata[currC][1] << ", " <<
                           idata[currC][2]);
            for (unsigned i = 0; i < idata.size(); i++) {
                // cannot reuse a, b, c
                if (i == curr.first || i == curr.second || i == currC) {
                    continue;
                }

                // compute directional distance
                dDist = dot(idata[i][0] - idata[curr.second][0],
                            idata[i][1] - idata[curr.second][1],
                            idata[i][2] - idata[curr.second][2],
                            vperp[0], vperp[1], vperp[2]);
                
                if (dDist > maxDist + EPS_LOC) {
                    currC = i;
                    maxDist = dDist;
                } else if (fabs(dDist) < EPS_LOC) {
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

            vcb = {idata[currC][0] - idata[curr.second][0],
                   idata[currC][1] - idata[curr.second][1],
                   idata[currC][2] - idata[curr.second][2]};
            // find new perpendicular vector
        } while (prevC != currC); // c not changed => all points on one side

        c = currC;

        R("FINAL C: " << idata[c][0] << ", " <<
                         idata[c][1] << ", " <<
                         idata[c][2]);

        // PROBLEM
        // this returns all points on faces, but I only need their convex hull
        // we need a way to find convex hull od points on plane, but in 3d
        
        // POSSIBLE SOLUTION
        // eliminate one non-zero coordinate
        
        onPlane.push_back(curr.first);
        onPlane.push_back(curr.second);

        R("PERPENDICULAR: " << vperp[0] << ", " << vperp[1] <<
                ", " << vperp[2]);
        Points2D planar; // 3d points on plane converted to 2d
        std::pair<int, int> coords; // id of coordinates to use
        if (fabs(vperp[0]) > EPS_LOC) {
            coords = {1, 2};
        } else if (fabs(vperp[1]) > EPS_LOC) {
            coords = {0, 2};
        } else {
            coords = {0, 1};
        }
        R("USING COORDS " << coords.first << ", " << coords.second);

        R("CONSTRUCTING FACE HULL OF");
        for (auto& i : onPlane) {
            R("  " << idata[i][0] << ", " << idata[i][1] << ", " << idata[i][2]
                   << " aka " <<
                      idata[i][coords.first] << ", " <<
                      idata[i][coords.second] << " (id: " << i << ")");
            planar.add({idata[i][coords.first], idata[i][coords.second]});
        }

        // find points on diameter of face in ccw order
        GrahamScan2D solver;
        std::vector<unsigned> faceID;
        solver.solveID(planar, faceID);
        unsigned fs = faceID.size();

        R("RESULTING FACE HULL")
        for (auto& i : faceID) {
            R("  " << idata[onPlane[i]][coords.first] << ", " <<
                      idata[onPlane[i]][coords.second] <<
                      " (id: " << onPlane[i] << ")");
        }

        bool flipped = 0;
        for (unsigned i = 0; i < fs; i++) {
            unsigned ex = onPlane[faceID[i]],
                     ey = onPlane[faceID[(i+1) % fs]];
            if (ex == curr.second && ey == curr.first) {
                R("EDGE ROTATION");
                flipped = 1;
            }
        }

        // add new found fresh to list, remove already found ones
        for (unsigned i = 0; i < fs; i++) {
            unsigned ex = onPlane[faceID[i]],
                     ey = onPlane[faceID[(i+1) % fs]];
            edge_t oe = {ex, ey};
            if (ex >= ey) {
                oe = {ey, ex};
            }
            if (flipped) {
                std::swap(ex, ey);
            }
            if (closed.find(oe) != closed.end()) {
                R("CLOSE EDGE: " << idata[ex][0] << ", " <<
                                    idata[ex][1] << ", " <<
                                    idata[ex][2] << " | " <<
                                    idata[ey][0] << ", " <<
                                    idata[ey][1] << ", " <<
                                    idata[ey][2] << " (id: " << ex <<
                                    ", " << ey << ")");
                fresh.erase({ey, ex});
                fresh.erase({ex, ey});
            } else {
                R("OPEN EDGE: " << idata[ey][0] << ", " <<
                                     idata[ey][1] << ", " <<
                                     idata[ey][2] << " | " <<
                                     idata[ex][0] << ", " <<
                                     idata[ex][1] << ", " <<
                                     idata[ex][2] << " (id: " << ex <<
                                    ", " << ey << ")");
                fresh.insert({ey, ex});
                closed.insert(oe);
            }
        }

        Points3D face;
        for (auto i : faceID) {
            face.add(idata[onPlane[i]]);
        }
        output.addFace(face);
    } while (!fresh.empty());

    // polyhedron with two faces means a plane => remove one of the faces
    if (output.getSize() == 2) {
        output.popFace();
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
            if (dif > EPS_LOC) {
                maxd = d;
                farcnt = 1;
                far = i;
            } else if (fabs(dif) < EPS_LOC) {
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
        if (dif < -EPS_LOC) {
            paired = i;
            maxd = d;
            minLen = iLen;
        } else if (fabs(dif) < EPS_LOC) {
            if (iLen < minLen - EPS_LOC) {
                paired = i;
                minLen = iLen;
            }
        }
    }

    return {far, paired};
}

}
