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
    const data_t& inputData = input.getData();
    std::pair<unsigned, unsigned> init = findInitial(inputData);
    R(init.first << ", " << init.second);
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

    do {
        // random direction vector
        far = 0;
        farcnt = 1;
        maxd = DBL_MIN;
        rndx = randomOne(); rndy = randomOne(); rndz = randomOne();
        R(rndx << ", " << rndy << ", " << rndz);
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
    for (unsigned i = 0; i < input.size(); i++) {
        if (i == far) {
            continue;
        }
        point_t vFarI = {input[far][0] - input[i][0],
                         input[far][1] - input[i][1],
                         input[far][2] - input[i][2]};
        double iLen = sqrt(vFarI[0]*vFarI[0] + vFarI[1]*vFarI[1]
                           + vFarI[2]*vFarI[2]);
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
