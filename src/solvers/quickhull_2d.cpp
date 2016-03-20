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
    return solveNaive(input, output);
}

void Quickhull2D::recNaive(point_t& a, point_t& b, data_t& plane)
{
    if (plane.size() == 0)
        return;

    // find point c furthest from ab
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

    std::pair<point_t, point_t> pivots = farthestPoints(inputData);
    point_t pivotLeft  = pivots.first,
            pivotRight = pivots.second;

    data_t topPlane, botPlane;
    for (auto& pt : inputData) {
        int o = orientation(pivotLeft[0],  pivotLeft[1],
                            pivotRight[0], pivotRight[1],
                            pt[0],         pt[1]);
        if (o == 1) {
            topPlane.push_back(pt);
        } else if (o == 2) {
            botPlane.push_back(pt);
        }
    }

    if (pivotLeft != pivotRight) {
        output.add(pivotLeft);
    }

    // recursive part
    recNaive(pivotLeft, pivotRight, topPlane);
    recNaive(pivotRight, pivotLeft, botPlane);

    output.add(pivotRight);

    return output;
}

Points2D& Quickhull2D::solveIterative(const Points2D& input, Points2D& output)
{
    if (input.getSize() <= 2) {
        output = input;
        return output;
    }
    const data_t& inputData = input.getData();

    std::pair<point_t, point_t> pivots;

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
            if (minX[1] - EPS > points[i][1]) {
                minX = points[i];
            }
        }

        delta = maxX[0] - points[i][0];
        if (delta < -EPS) {
            maxX = points[i];
        } else if (fabs(delta) < EPS) {
            if (maxX[1] + EPS < points[i][1]) {
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

    std::pair<point_t, point_t> furthest = {candidates[0], candidates[1]};
    double maxDist = 0.0, currDist;
    for (int i = 0; i < 4; i++) {
        for (int j = i + 1; j < 4; j++) {
            currDist = dist({candidates[i][0], candidates[i][1]},
                            {candidates[j][0], candidates[j][1]});
            if (currDist > maxDist) {
                furthest = {candidates[i], candidates[j]};
                maxDist = currDist;
            }
        }
    }

    if (furthest.first[0] < furthest.second[0]) {
        std::swap(furthest.first, furthest.second);
    }
    return furthest;
}

}
