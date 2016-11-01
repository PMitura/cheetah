#include "solvers/quickhull_2d.h"

namespace ch
{

Quickhull2D::Quickhull2D()
    :globIn_(NULL), globOut_(NULL)
{
    name_ = "Quickhull";
    EPS_LOC = 1e-6;
    variant_ = FORWARD;
    parallelThreshold_ = 1000;
}

Quickhull2D::Quickhull2D(Variant v)
    :globIn_(NULL), globOut_(NULL)
{
    name_ = "Quickhull";
    EPS_LOC = 1e-6;
    variant_ = v;
    parallelThreshold_ = 1000;
}

Quickhull2D::Quickhull2D(Variant v, int threshold)
    :globIn_(NULL), globOut_(NULL)
{
    name_ = "Quickhull";
    EPS_LOC = 1e-6;
    variant_ = v;
    parallelThreshold_ = threshold;
}

Points2D& Quickhull2D::solve(const Points2D& input, Points2D& output)
{
    // temp lower global eps
    // EPS = 1e-6;

    switch (variant_) {
        case NAIVE:
            solveSequential(input, output);
            break;
        case PRECOMP:
            solvePrecomp(input, output);
            break;
        case FORWARD:
            solveForwarded(input, output);
            break;
        case PARA:
            solveParallel(input, output);
            break;
    }

    // EPS = 1e-12;
    return output;

    // discarded variants
    // return solveNaive(input, output);
    // return solveIterative(input, output);
}

void Quickhull2D::recNaive(point_t& a, point_t& b, data_t& plane)
{
    if (plane.size() == 0)
        return;

    // find point c farthest from ab
    point_t c = planeFarthestDist(a, b, plane);

    data_t acPlane, cbPlane;

    for (auto& pt : plane) {
        if (orientHiEPS(a[0], a[1], c[0], c[1], pt[0], pt[1]) == 2) {
            acPlane.push_back(pt);
        } else if (orientHiEPS(c[0], c[1], b[0], b[1], pt[0], pt[1]) == 2) {
            cbPlane.push_back(pt);
        }
    }

    recNaive(a, c, acPlane);
    globOut_ -> add(c);
    recNaive(c, b, cbPlane);
}

Points2D& Quickhull2D::solveNaive(const Points2D& input, Points2D& output)
{
    if (input.getSize() <= 2) {
        output = input;
        return output;
    }

    globOut_ = &output;
    const data_t& inputData = input.getData();

    // alt std::pair<point_t, point_t> pivots = minMaxX(inputData);
    std::pair<point_t, point_t> pivots = farthestPoints(inputData);
    point_t pivotLeft  = pivots.first,
            pivotRight = pivots.second;

    data_t topPlane, botPlane;
    divideToPlanes(inputData, pivotLeft, pivotRight, topPlane, botPlane);

    // recursive part
    output.add(pivotLeft);
    recNaive(pivotLeft, pivotRight, topPlane);
    output.add(pivotRight);
    recNaive(pivotRight, pivotLeft, botPlane);

    return output;
}


void Quickhull2D::recSplit(point_t& a, point_t& b, point_t& c,
                           data_t& plane, bool upper)
{
    if (plane.size() == 0)
        return;

    data_t acPlane, cbPlane;
    double acMax = -1, cbMax = -1;
    point_t acFar, cbFar;

    double aco, cbo;
    if (upper) {
        for (auto& pt : plane) {
            if (pt[0] < c[0] - EPS_LOC) {
                aco = cross(a[0], a[1], c[0], c[1], pt[0], pt[1]);
                if (aco < -EPS_LOC) {
                    acPlane.push_back(pt);
                    if (fabs(aco) > acMax) {
                        acFar = pt;
                        acMax = fabs(aco);
                    }
                    continue;
                }
            }

            if (pt[0] > c[0] + EPS_LOC) {
                cbo = cross(c[0], c[1], b[0], b[1], pt[0], pt[1]);
                if (cbo < -EPS_LOC) {
                    cbPlane.push_back(pt);
                    if (fabs(cbo) > cbMax) {
                        cbFar = pt;
                        cbMax = fabs(cbo);
                    }
                }
            }
        }
    } else {
        for (auto& pt : plane) {
            if (pt[0] > c[0] + EPS_LOC) {
                aco = cross(a[0], a[1], c[0], c[1], pt[0], pt[1]);
                if (aco < -EPS_LOC) {
                    acPlane.push_back(pt);
                    if (fabs(aco) > acMax) {
                        acFar = pt;
                        acMax = fabs(aco);
                    }
                    continue;
                }
            }

            if (pt[0] < c[0] - EPS_LOC) {
                cbo = cross(c[0], c[1], b[0], b[1], pt[0], pt[1]);
                if (cbo < -EPS_LOC) {
                    cbPlane.push_back(pt);
                    if (fabs(cbo) > cbMax) {
                        cbFar = pt;
                        cbMax = fabs(cbo);
                    }
                }
            }
        }
    }


    recSplit(a, c, acFar, acPlane, upper);
    globOut_ -> add(c);
    recSplit(c, b, cbFar, cbPlane, upper);
}

void Quickhull2D::recSequential(point_t& a, point_t& b, point_t& c,
                             data_t& plane)
{
    if (plane.size() == 0)
        return;

    data_t acPlane, cbPlane;
    double acMax = -1, cbMax = -1;
    point_t acFar, cbFar;

    double aco, cbo;
    for (auto& pt : plane) {
        aco = cross(a[0], a[1], c[0], c[1], pt[0], pt[1]);
        if (aco < -EPS_LOC) {
            acPlane.push_back(pt);
            if (fabs(aco) > acMax) {
                acFar = pt;
                acMax = fabs(aco);
            }
            continue;
        }

        cbo = cross(c[0], c[1], b[0], b[1], pt[0], pt[1]);
        if (cbo < -EPS_LOC) {
            cbPlane.push_back(pt);
            if (fabs(cbo) > cbMax) {
                cbFar = pt;
                cbMax = fabs(cbo);
            }
        }
    }

    recSequential(a, c, acFar, acPlane);
    globOut_ -> add(c);
    recSequential(c, b, cbFar, cbPlane);
}

Points2D& Quickhull2D::solveSequential(const Points2D& input, Points2D& output)
{
    if (input.getSize() <= 2) {
        output = input;
        return output;
    }

    globOut_ = &output;
    const data_t& inputData = input.getData();

    std::pair<point_t, point_t> pivots = minMaxX(inputData);
    // std::pair<point_t, point_t> pivots = farthestPoints(inputData);
    point_t pivotLeft  = pivots.first,
            pivotRight = pivots.second;

    data_t topPlane, botPlane;
    // future farthest points
    double topMax = -1, botMax = -1;
    point_t topFar, botFar;

    /*
    // precompute cross product
    double alpha = pivotLeft[1] - pivotRight[1],
           beta  = pivotLeft[0] - pivotRight[0],
           gamma = beta*pivotLeft[1] - alpha*pivotLeft[0];
           */

    // extended divide to planes
    for (auto& pt : inputData) {
        double o = cross(pivotLeft[0],  pivotLeft[1],
                         pivotRight[0], pivotRight[1],
                         pt[0],         pt[1]);
        if (o < -EPS_LOC) {
            topPlane.push_back(pt);
            if (fabs(o) > topMax) {
                topFar = pt;
                topMax = fabs(o);
            }
        } else if (o > EPS_LOC) {
            botPlane.push_back(pt);
            if (fabs(o) > botMax) {
                botFar = pt;
                botMax = fabs(o);
            }
        }
    }

    // recursive part
    output.add(pivotLeft);
    recSequential(pivotLeft, pivotRight, topFar, topPlane);
    output.add(pivotRight);
    recSequential(pivotRight, pivotLeft, botFar, botPlane);

    return output;
}

void Quickhull2D::recPrecomp(point_t& a, point_t& b, point_t& c,
                             data_t& plane)
{
    if (plane.size() == 0)
        return;

    data_t acPlane, cbPlane;
    double acMax = -1, cbMax = -1;
    point_t acFar, cbFar;

    double aco, cbo;
    // precompute cross
    double alphaAC = a[1] - c[1],
           betaAC  = a[0] - c[0],
           gammaAC = betaAC*c[1] - alphaAC*c[0],
           alphaCB = c[1] - b[1],
           betaCB  = c[0] - b[0],
           gammaCB = betaCB*b[1] - alphaCB*b[0];
    for (auto& pt : plane) {
        aco = partCross(pt[0], pt[1], alphaAC, betaAC, gammaAC);
        if (aco < -EPS_LOC) {
            acPlane.push_back(pt);
            if (fabs(aco) > acMax) {
                acFar = pt;
                acMax = fabs(aco);
            }
            continue;
        }

        cbo = partCross(pt[0], pt[1], alphaCB, betaCB, gammaCB);
        if (cbo < -EPS_LOC) {
            cbPlane.push_back(pt);
            if (fabs(cbo) > cbMax) {
                cbFar = pt;
                cbMax = fabs(cbo);
            }
        }
    }

    recPrecomp(a, c, acFar, acPlane);
    globOut_ -> add(c);
    recPrecomp(c, b, cbFar, cbPlane);
}

Points2D& Quickhull2D::solvePrecomp(const Points2D& input, Points2D& output)
{
    if (input.getSize() <= 2) {
        output = input;
        return output;
    }

    globOut_ = &output;
    const data_t& inputData = input.getData();

    std::pair<point_t, point_t> pivots = minMaxX(inputData);
    // std::pair<point_t, point_t> pivots = farthestPoints(inputData);
    point_t pivotLeft  = pivots.first,
            pivotRight = pivots.second;

    data_t topPlane, botPlane;
    // future farthest points
    double topMax = -1, botMax = -1;
    point_t topFar, botFar;

    // precompute cross product
    double alpha = pivotLeft[1] - pivotRight[1],
           beta  = pivotLeft[0] - pivotRight[0],
           gamma = beta*pivotLeft[1] - alpha*pivotLeft[0];

    // extended divide to planes
    for (auto& pt : inputData) {
        double o = partCross(pt[0], pt[1], alpha, beta, gamma);
        if (o < -EPS_LOC) {
            topPlane.push_back(pt);
            if (fabs(o) > topMax) {
                topFar = pt;
                topMax = fabs(o);
            }
        } else if (o > EPS_LOC) {
            botPlane.push_back(pt);
            if (fabs(o) > botMax) {
                botFar = pt;
                botMax = fabs(o);
            }
        }
    }

    // recursive part
    output.add(pivotLeft);
    recPrecomp(pivotLeft, pivotRight, topFar, topPlane);
    output.add(pivotRight);
    recPrecomp(pivotRight, pivotLeft, botFar, botPlane);

    return output;
}

void Quickhull2D::recForwarded(const point_t& a, const point_t& b, const point_t& c,
                              std::vector<unsigned> plane, unsigned planeSize)
{
    if (planeSize == 0) {
        return;
    }

    std::vector<unsigned> acPlane, cbPlane;
    double acMax = -1, cbMax = -1;
    unsigned acFar = 0, cbFar = 0;

    /*
    // precompute angle
    double alphaAC = a[1] - c[1],
           betaAC  = a[0] - c[0],
           gammaAC = betaAC*c[1] - alphaAC*c[0],
           alphaCB = c[1] - b[1],
           betaCB  = c[0] - b[0],
           gammaCB = betaCB*b[1] - alphaCB*b[0];
    double * acPre = new double[plane.size()],
           * cbPre = new double[plane.size()];
    for (unsigned i = 0; i < planeSize; i++) {
        acPre[i] = partCross((*globIn_)[pt][0], (*globIn_)[pt][1],
                             alphaAC, betaAC, gammaAC);
        cbPre[i] = partCross((*globIn_)[pt][0], (*globIn_)[pt][1],
                             alphaCB, betaCB, gammaCB);
    }
    */

    for (unsigned i = 0; i < planeSize; i++) {
        int pt = plane[i];
        double aco = cross(a[0], a[1], c[0], c[1],
                        (*globIn_)[pt][0], (*globIn_)[pt][1]);
        if (aco > EPS_LOC) {
            acPlane.push_back(pt);
            double fac = fabs(aco);
            if (fac > acMax) {
                acFar = pt;
                acMax = fac;
            }
            continue;
        }

        double cbo = cross(c[0], c[1], b[0], b[1],
                        (*globIn_)[pt][0], (*globIn_)[pt][1]);
        if (cbo > EPS_LOC) {
            cbPlane.push_back(pt);
            double fcb = fabs(cbo);
            if (fcb > cbMax) {
                cbFar = pt;
                cbMax = fcb;
            }
        }
    }

    recForwarded(a, c, (*globIn_)[acFar], acPlane, acPlane.size());
    globOut_ -> add(c);
    recForwarded(c, b, (*globIn_)[cbFar], cbPlane, cbPlane.size());
}

Points2D& Quickhull2D::solveForwarded(const Points2D& input, Points2D& output)
{
    if (input.getSize() <= 2) {
        output = input;
        return output;
    }

    globOut_ = &output;
    const data_t& inputData = input.getData();
    globIn_ = &(input.getData());

    std::pair<point_t, point_t> pivots = minMaxX(inputData);
    // std::pair<point_t, point_t> pivots = farthestPoints(inputData);
    point_t pivotLeft  = pivots.first,
            pivotRight = pivots.second;

    std::vector<unsigned> topPlane, botPlane;
    topPlane.resize(input.getSize());
    botPlane.resize(input.getSize());
    unsigned topPtr = 0, botPtr = 0;
    // future farthest points
    double topMax = -1, botMax = -1;
    unsigned topFar = 0, botFar = 0;

    /*
    // precompute cross product
    double alpha = pivotLeft[1] - pivotRight[1],
           beta  = pivotLeft[0] - pivotRight[0],
           gamma = beta*pivotRight[1] - alpha*pivotRight[0];
    // precomputed cross has lower precision
    double oldEPS_LOC = EPS_LOC;
    EPS_LOC = 1e-6; 
    */

    // extended divide to planes
    for (unsigned i = 0; i < inputData.size(); i++) {
        // double o = partCross(inputData[i][0], inputData[i][1],
        double o = cross(pivotLeft[0],    pivotLeft[1],
                         pivotRight[0],   pivotRight[1],
                         inputData[i][0], inputData[i][1]);
        double oa = fabs(o);
        if (o < -EPS_LOC) {
            topPlane[topPtr++] = i;
            if (oa > topMax) {
                topFar = i;
                topMax = oa;
            }
        } else if (o > EPS_LOC) {
            botPlane[botPtr++] = i;
            if (oa > botMax) {
                botFar = i;
                botMax = oa;
            }
        }
    }

    output.add(pivotRight);
    recForwarded(pivotRight, pivotLeft, inputData[topFar], topPlane,
                topPtr);
    output.add(pivotLeft);
    recForwarded(pivotLeft, pivotRight, inputData[botFar], botPlane,
                botPtr);

    return output;
}

void Quickhull2D::recParallel(const point_t& a, const point_t& b, const point_t& c,
                              std::vector<unsigned> plane, unsigned planeSize,
                              std::list<point_t>& onHull)
{
    if (planeSize == 0) {
        return;
    }

    std::vector<unsigned> acPlane, cbPlane;
    double acMax = -1, cbMax = -1;
    unsigned acFar = 0, cbFar = 0;

    
    /*
    // precompute angle
    double alphaAC = a[1] - c[1],
           betaAC  = a[0] - c[0],
           gammaAC = betaAC*c[1] - alphaAC*c[0],
           alphaCB = c[1] - b[1],
           betaCB  = c[0] - b[0],
           gammaCB = betaCB*b[1] - alphaCB*b[0];
    double * acPre = new double[plane.size()],
           * cbPre = new double[plane.size()];
    for (unsigned i = 0; i < planeSize; i++) {
        acPre[i] = partCross((*globIn_)[pt][0], (*globIn_)[pt][1],
                             alphaAC, betaAC, gammaAC);
        cbPre[i] = partCross((*globIn_)[pt][0], (*globIn_)[pt][1],
                             alphaCB, betaCB, gammaCB);
    }
    */

    for (unsigned i = 0; i < planeSize; i++) {
        int pt = plane[i];
        double aco = cross(a[0], a[1], c[0], c[1],
                        (*globIn_)[pt][0], (*globIn_)[pt][1]);
        if (aco > EPS_LOC) {
            acPlane.push_back(pt);
            double fac = fabs(aco);
            if (fac > acMax) {
                acFar = pt;
                acMax = fac;
            }
            continue;
        }

        double cbo = cross(c[0], c[1], b[0], b[1],
                        (*globIn_)[pt][0], (*globIn_)[pt][1]);
        if (cbo > EPS_LOC) {
            cbPlane.push_back(pt);
            double fcb = fabs(cbo);
            if (fcb > cbMax) {
                cbFar = pt;
                cbMax = fcb;
            }
        }
    }

    // delete[] acPre;
    // delete[] cbPre;

    std::list<point_t> acList, cbList;

    if (acPlane.size() > parallelThreshold_ &&
        cbPlane.size() > parallelThreshold_) {
        // R("PARASPLIT " << acPlane.size() << " - " << cbPlane.size());
#pragma omp parallel
        {
#pragma omp sections
            {
#pragma omp section
                {
                recParallel(a, c, (*globIn_)[acFar], acPlane, acPlane.size(),
                            acList);
                }
#pragma omp section
                {
                recParallel(c, b, (*globIn_)[cbFar], cbPlane, cbPlane.size(),
                            cbList);
                }
            }
        }
    } else {
        recParallel(a, c, (*globIn_)[acFar], acPlane, acPlane.size(), acList);
        recParallel(c, b, (*globIn_)[cbFar], cbPlane, cbPlane.size(), cbList);
    }


    // O(1) append to onHull
    onHull.splice(onHull.end(), acList);
    onHull.push_back(c);
    onHull.splice(onHull.end(), cbList);
}

Points2D& Quickhull2D::solveParallel(const Points2D& input, Points2D& output)
{
    if (input.getSize() <= 2) {
        output = input;
        return output;
    }

    const data_t& inputData = input.getData();
    globIn_ = &(input.getData());
    parallelThreshold_ = 100;

    std::pair<point_t, point_t> pivots = minMaxX(inputData);
    // std::pair<point_t, point_t> pivots = farthestPoints(inputData);
    point_t pivotLeft  = pivots.first,
            pivotRight = pivots.second;


    std::list<point_t> topList, botList;
    std::vector<unsigned> topPlane, botPlane;
    topPlane.resize(input.getSize());
    botPlane.resize(input.getSize());
    unsigned topPtr = 0, botPtr = 0;
    // future farthest points
    double topMax = -1, botMax = -1;
    unsigned topFar = 0, botFar = 0;

    /*
    // precompute cross product
    double alpha = pivotLeft[1] - pivotRight[1],
           beta  = pivotLeft[0] - pivotRight[0],
           gamma = beta*pivotRight[1] - alpha*pivotRight[0];
    // precomputed cross has lower precision
    double oldEPS_LOC = EPS_LOC;
    EPS_LOC = 1e-6; 
    */

    // extended divide to planes

    std::vector<double> crosses;
    crosses.resize(input.getSize());
#pragma omp parallel for default(shared) schedule(static)
    for (unsigned i = 0; i < inputData.size(); i++) {
        crosses[i] = cross(pivotLeft[0],    pivotLeft[1],
                           pivotRight[0],   pivotRight[1],
                           inputData[i][0], inputData[i][1]);
    }
    for (unsigned i = 0; i < inputData.size(); i++) {
        double o = crosses[i];
        double oa = fabs(o);
        if (o < -EPS_LOC) {
            topPlane[topPtr++] = i;
            if (oa > topMax) {
                topFar = i;
                topMax = oa;
            }
        } else if (o > EPS_LOC) {
            botPlane[botPtr++] = i;
            if (oa > botMax) {
                botFar = i;
                botMax = oa;
            }
        }
    }

    /*
    recParallel(pivotLeft, pivotRight, inputData[topFar], topPlane,
                topPtr, topList);
    recParallel(pivotRight, pivotLeft, inputData[botFar], botPlane,
                botPtr, botList);
                */

    omp_set_dynamic(1);
    omp_set_nested(5);
#pragma omp parallel
    {
#pragma omp sections
        {
#pragma omp section
            recParallel(pivotRight, pivotLeft, inputData[topFar], topPlane,
                        topPtr, topList);
#pragma omp section
            recParallel(pivotLeft, pivotRight, inputData[botFar], botPlane,
                        botPtr, botList);
        }
    }

    output.add(pivotRight);
    for (auto pt : topList) {
        output.add(pt);
    }
    output.add(pivotLeft);
    for (auto pt : botList) {
        output.add(pt);
    }

    // EPS_LOC = oldEPS_LOC; // return EPS_LOC back to previous state

    return output;
}

Points2D& Quickhull2D::solveIterative(const Points2D& input, Points2D& output)
{
    if (input.getSize() <= 2) {
        output = input;
        return output;
    }
    const data_t& inputData = input.getData();

    struct Face {
        point_t a, b;
        data_t see;
        bool save;
    };

    // alt std::pair<point_t, point_t> pivots = minMaxX(inputData);
    std::pair<point_t, point_t> pivots = farthestPoints(inputData);
    point_t pivotLeft  = pivots.first,
            pivotRight = pivots.second;

    Face * topFace = new Face, * botFace = new Face;
    topFace -> a = pivotLeft;  topFace -> b = pivotRight; topFace -> save = 0;
    botFace -> a = pivotRight; botFace -> b = pivotLeft;  botFace -> save = 0;
    divideToPlanes(inputData, pivotLeft, pivotRight,
                   topFace -> see, botFace -> see);

    std::vector<Face*> bases;
    bases.push_back(topFace);
    bases.push_back(botFace);

    for (auto base : bases) {
        output.add(base -> a);
        std::stack<Face*> faces;
        faces.push(base);

        while (!faces.empty()) {
            Face * curr = faces.top();
            faces.pop();
            if (curr -> save) {
                output.add(curr -> a);
                delete curr;
                continue;
            }
            if (curr -> see.size() == 0) {
                delete curr;
                continue;
            }

            point_t a = curr -> a, b = curr -> b,
                    c = planeFarthestDist(a, b, curr -> see);
            Face * ac = new Face, * cb = new Face, * adder = new Face;
            ac -> a = a; ac -> b = c; ac -> save = 0;
            cb -> a = c; cb -> b = b; cb -> save = 0;
            for (auto& pt : curr -> see) {
                if (orientHiEPS(a[0], a[1], c[0], c[1], pt[0], pt[1]) == 1) {
                    ac -> see.push_back(pt);
                } else if (orientHiEPS(c[0], c[1], b[0], b[1],
                                       pt[0], pt[1]) == 1) {
                    cb -> see.push_back(pt);
                }
            }
            adder -> a = c;
            adder -> save = 1;

            faces.push(ac);
            faces.push(adder);
            faces.push(cb);

            delete curr;
        }
    }

    return output;
}

std::pair<point_t, point_t> Quickhull2D::minMaxX(const data_t& points)
{
    double delta;
    point_t minX = points[0], maxX = points[0];
    for (unsigned i = 1; i < points.size(); i++) {
        delta = minX[0] - points[i][0];
        if (delta > EPS_LOC) {
            minX = points[i];
        } else if (fabs(delta) < EPS_LOC) {
            if (minX[1] + EPS_LOC < points[i][1]) {
                minX = points[i];
            }
        }

        delta = maxX[0] - points[i][0];
        if (delta < -EPS_LOC) {
            maxX = points[i];
        } else if (fabs(delta) < EPS_LOC) {
            if (maxX[1] - EPS_LOC > points[i][1]) {
                maxX = points[i];
            }
        }
    }

    return {minX, maxX};
}

std::pair<point_t, point_t> Quickhull2D::farthestPoints(const data_t& points)
{
    point_t minX = points[0], maxX = points[0],
            minY = points[0], maxY = points[0];
    for (unsigned i = 1; i < points.size(); i++) {
        if (points[i][0] < minX[0]) {
            minX = points[i];
        }
        if (points[i][0] > maxX[0]) {
            maxX = points[i];
        }
        if (points[i][1] < minY[1]) {
            minY = points[i];
        }
        if (points[i][1] > maxY[1]) {
            maxY = points[i];
        }
    }
    std::vector<point_t> candidates;
    candidates.push_back(minY);
    candidates.push_back(maxY);
    candidates.push_back(minX);
    candidates.push_back(minY);

    std::pair<point_t, point_t> farthest = {candidates[0], candidates[1]};
    double maxDist = 0.0, currDist;
    for (int i = 0; i < 4; i++) {
        for (int j = i + 1; j < 4; j++) {
            currDist = dist({candidates[i][0], candidates[i][1]},
                            {candidates[j][0], candidates[j][1]});
            if (currDist > maxDist) {
                farthest = {candidates[i], candidates[j]};
                maxDist = currDist;
            }
        }
    }

    if (farthest.first[0] < farthest.second[0]) {
        std::swap(farthest.first, farthest.second);
    }
    return farthest;
}

point_t Quickhull2D::planeFarthestCross(point_t& a, point_t& b,
                                        const data_t& plane)
{
    point_t c = plane[0];
    double maxCross = cross(a[0], a[1], b[0], b[1], c[0], c[1]),
           currCross;
    for (auto& pt : plane) {
        currCross = cross(a[0], a[1], b[0], b[1], pt[0], pt[1]);
        if (maxCross - currCross < -EPS_LOC) {
            maxCross = currCross;
            c = pt;
        }
    }
    return c;
}

point_t Quickhull2D::planeFarthestDist(point_t& a, point_t& b,
                                        const data_t& plane)
{
    point_t c = plane[0];
    double maxDist = distToLine({a[0], a[1]}, {b[0], b[1]}, {c[0], c[1]}),
           currDist;
    for (auto& pt : plane) {
        currDist = distToLine({a[0], a[1]}, {b[0], b[1]}, {pt[0], pt[1]});
        if (maxDist - currDist < -EPS_LOC) {
            maxDist = currDist;
            c = pt;
        }
    }
    return c;
}

void Quickhull2D::divideToPlanes(const data_t& input,
                                 point_t& pivotLeft, point_t& pivotRight,
                                 data_t& topPlane, data_t& botPlane)
{
    for (auto& pt : input) {
        int o = orientHiEPS(pivotLeft[0],  pivotLeft[1],
                            pivotRight[0], pivotRight[1],
                            pt[0],         pt[1]);
        if (o == 2) {
            topPlane.push_back(pt);
        } else if (o == 1) {
            botPlane.push_back(pt);
        }
    }
}

void Quickhull2D::divideToPlanesPara(const data_t& input,
                                     point_t& pivotLeft, point_t& pivotRight,
                                     data_t& topPlane, data_t& botPlane)
{
    int * med = new int[input.size()];

#pragma omp parallel for default(shared) schedule(static)
    for (unsigned i = 0; i < input.size(); i++) {
        med[i] = orientHiEPS(pivotLeft[0],  pivotLeft[1],
                             pivotRight[0], pivotRight[1],
                             input[i][0],   input[i][1]);
    }

    for (unsigned i = 0; i < input.size(); i++) {
        if (med[i] == 2) {
            topPlane.push_back(input[i]);
        } else if (med[i] == 1) {
            botPlane.push_back(input[i]);
        }
    }

    delete[] med;
}

}
