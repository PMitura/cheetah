#include "solvers/quickhull_2d.h"

namespace ch
{

Quickhull2D::Quickhull2D()
    :globOut_(NULL)
{
    name_ = "Quickhull";
}

Points2D& Quickhull2D::solve(const Points2D& input, Points2D& output)
{
    // return solveParallel(input, output);
    // return solveNaive(input, output);
    return solveOptimal(input, output);
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
        if (orientation(a[0], a[1], c[0], c[1], pt[0], pt[1]) == 1) {
            acPlane.push_back(pt);
        } else if (orientation(c[0], c[1], b[0], b[1], pt[0], pt[1]) == 1) {
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

void Quickhull2D::recParallel(point_t a, point_t b, data_t& plane,
                              std::list<point_t>& onHull)
{
    if (plane.size() == 0) {
        return;
    }

    // find point c farthest from ab
    point_t c = planeFarthestDist(a, b, plane);

    data_t acPlane, cbPlane;

    for (auto& pt : plane) {
        if (orientation(a[0], a[1], c[0], c[1], pt[0], pt[1]) == 1) {
            acPlane.push_back(pt);
        } else if (orientation(c[0], c[1], b[0], b[1], pt[0], pt[1]) == 1) {
            cbPlane.push_back(pt);
        }
    }
    std::list<point_t> acList, cbList;

    bool cond = (acPlane.size() > plane.size() / 10) && plane.size() > 1000;
    /*
    if (cond) {
        R("split @" << acPlane.size())
    }
    */

#pragma omp task shared(acPlane, acList) if (cond)
        recParallel(a, c, acPlane, acList);
    recParallel(c, b, cbPlane, cbList);
#pragma omp taskwait

    // O(1) append to onHull
    onHull.splice(onHull.end(), acList);
    onHull.push_back(c);
    onHull.splice(onHull.end(), cbList);
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
            if (pt[0] < c[0] - EPS) {
                aco = cross(a[0], a[1], c[0], c[1], pt[0], pt[1]);
                if (aco > EPS) {
                    acPlane.push_back(pt);
                    if (fabs(aco) > acMax) {
                        acFar = pt;
                        acMax = fabs(aco);
                    }
                    continue;
                }
            }

            if (pt[0] > c[0] + EPS) {
                cbo = cross(c[0], c[1], b[0], b[1], pt[0], pt[1]);
                if (cbo > EPS) {
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
            if (pt[0] > c[0] + EPS) {
                aco = cross(a[0], a[1], c[0], c[1], pt[0], pt[1]);
                if (aco > EPS) {
                    acPlane.push_back(pt);
                    if (fabs(aco) > acMax) {
                        acFar = pt;
                        acMax = fabs(aco);
                    }
                    continue;
                }
            }

            if (pt[0] < c[0] - EPS) {
                cbo = cross(c[0], c[1], b[0], b[1], pt[0], pt[1]);
                if (cbo > EPS) {
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

void Quickhull2D::recOptimal(point_t& a, point_t& b, point_t& c,
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
        if (aco > EPS) {
            acPlane.push_back(pt);
            if (fabs(aco) > acMax) {
                acFar = pt;
                acMax = fabs(aco);
            }
            continue;
        }

        cbo = cross(c[0], c[1], b[0], b[1], pt[0], pt[1]);
        if (cbo > EPS) {
            cbPlane.push_back(pt);
            if (fabs(cbo) > cbMax) {
                cbFar = pt;
                cbMax = fabs(cbo);
            }
        }
    }

    recOptimal(a, c, acFar, acPlane);
    globOut_ -> add(c);
    recOptimal(c, b, cbFar, cbPlane);
}

Points2D& Quickhull2D::solveOptimal(const Points2D& input, Points2D& output)
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
        // double o = partCross(pt[0], pt[1], alpha, beta, gamma);
        if (o > EPS) {
            topPlane.push_back(pt);
            if (fabs(o) > topMax) {
                topFar = pt;
                topMax = fabs(o);
            }
        } else if (o < -EPS) {
            botPlane.push_back(pt);
            if (fabs(o) > botMax) {
                botFar = pt;
                botMax = fabs(o);
            }
        }
    }

    // recursive part
    output.add(pivotLeft);
    recOptimal(pivotLeft, pivotRight, topFar, topPlane);
    output.add(pivotRight);
    recOptimal(pivotRight, pivotLeft, botFar, botPlane);

    return output;
}

Points2D& Quickhull2D::solvePreprocessed(const Points2D& input,
                                         Points2D& output)
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
    // future farthest points
    double topMax = -1, botMax = -1;
    point_t topFar, botFar;

    // extended divide to planes
    for (auto& pt : inputData) {
        double o = cross(pivotLeft[0],  pivotLeft[1],
                         pivotRight[0], pivotRight[1],
                         pt[0],         pt[1]);
        // double o = partCross(pt[0], pt[1], alpha, beta, gamma);
        if (o > EPS) {
            topPlane.push_back(pt);
            if (fabs(o) > topMax) {
                topFar = pt;
                topMax = fabs(o);
            }
        } else if (o < -EPS) {
            botPlane.push_back(pt);
            if (fabs(o) > botMax) {
                botFar = pt;
                botMax = fabs(o);
            }
        }
    }

    // recursive part
    output.add(pivotLeft);
    // recOptimal(pivotLeft, pivotRight, topFar, topPlane);
    recSplit(pivotLeft, pivotRight, topFar, topPlane, true);
    output.add(pivotRight);
    // recOptimal(pivotRight, pivotLeft, botFar, botPlane);
    recSplit(pivotRight, pivotLeft, botFar, botPlane, false);

    return output;
}

Points2D& Quickhull2D::solveParallel(const Points2D& input, Points2D& output)
{
    if (input.getSize() <= 2) {
        output = input;
        return output;
    }

    const data_t& inputData = input.getData();

    // alt std::pair<point_t, point_t> pivots = minMaxX(inputData);
    std::pair<point_t, point_t> pivots = farthestPoints(inputData);
    point_t pivotLeft  = pivots.first,
            pivotRight = pivots.second;

    data_t topPlane, botPlane;
    divideToPlanes(inputData, pivotLeft, pivotRight, topPlane, botPlane);

    std::list<point_t> topList, botList;

#pragma omp parallel
    {
#pragma omp single
        {
#pragma omp task shared(topList, topPlane) if (topPlane.size() > input.getSize() / 3)
                recParallel(pivotLeft, pivotRight, topPlane, topList);
            recParallel(pivotRight, pivotLeft, botPlane, botList);
#pragma omp taskwait
        }
    }

    output.add(pivotLeft);
    for (auto pt : topList) {
        output.add(pt);
    }
    output.add(pivotRight);
    for (auto pt : botList) {
        output.add(pt);
    }

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
                if (orientation(a[0], a[1], c[0], c[1], pt[0], pt[1]) == 1) {
                    ac -> see.push_back(pt);
                } else if (orientation(c[0], c[1], b[0], b[1],
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
        if (delta > EPS) {
            minX = points[i];
        } else if (fabs(delta) < EPS) {
            if (minX[1] + EPS < points[i][1]) {
                minX = points[i];
            }
        }

        delta = maxX[0] - points[i][0];
        if (delta < -EPS) {
            maxX = points[i];
        } else if (fabs(delta) < EPS) {
            if (maxX[1] - EPS > points[i][1]) {
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
        if (maxCross - currCross < -EPS) {
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
        if (maxDist - currDist < -EPS) {
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
        int o = orientation(pivotLeft[0],  pivotLeft[1],
                            pivotRight[0], pivotRight[1],
                            pt[0],         pt[1]);
        if (o == 1) {
            topPlane.push_back(pt);
        } else if (o == 2) {
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
        med[i] = orientation(pivotLeft[0],  pivotLeft[1],
                             pivotRight[0], pivotRight[1],
                             input[i][0],   input[i][1]);
    }

    for (unsigned i = 0; i < input.size(); i++) {
        if (med[i] == 1) {
            topPlane.push_back(input[i]);
        } else if (med[i] == 2) {
            botPlane.push_back(input[i]);
        }
    }

    delete[] med;
}

}
