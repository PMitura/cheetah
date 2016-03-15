#include "solvers/quickhull_2d.h"

namespace ch
{

Quickhull2D::Quickhull2D()
    :globOut_(NULL)
{}

Points2D& Quickhull2D::solve(const Points2D& input, Points2D& output)
{
    return solveNaive(input, output);
}

void Quickhull2D::recNaive(point_t& a, point_t& b, data_t* plane)
{
    if (plane -> size() == 0)
        return;

    // find point c furthest from ab
    point_t c = (*plane)[0];
    double maxDist = distToLine({a[0], a[1]}, {b[0], b[1]}, {c[0], c[1]}),
           currDist;
    for (auto& pt : (*plane)) {
        currDist = distToLine({a[0], a[1]}, {b[0], b[1]}, {pt[0], pt[1]});
        if (maxDist - currDist < -EPS) {
            maxDist = currDist;
            c = pt;
        }
    }

    data_t* acPlane = new data_t, * cbPlane = new data_t;

    for (auto& pt : (*plane)) {
        if (orientation(a[0], a[1], c[0], c[1], pt[0], pt[1]) == 1) {
            acPlane -> push_back(pt);
        } else if (orientation(c[0], c[1], b[0], b[1], pt[0], pt[1]) == 1) {
            cbPlane -> push_back(pt);
        }
    }
    /*
    D(acPlane -> size());
    D(cbPlane -> size());
    */

    globOut_ -> add(c);
    recNaive(a, c, acPlane);
    recNaive(c, b, cbPlane);

    delete acPlane;
    delete cbPlane;
}

Points2D& Quickhull2D::solveNaive(const Points2D& input, Points2D& output)
{
    if (input.getSize() <= 2) {
        output = input;
        return output;
    }

    globOut_ = &output;
    const data_t& inputData = input.getData();

    // find min and max X, min/max Y in case of tie
    point_t minX = inputData[0], maxX = inputData[0];
    double delta;
    for (unsigned i = 1; i < inputData.size(); i++) {
        delta = minX[0] - inputData[i][0];
        if (delta > EPS) {
            minX = inputData[i];
        } else if (fabs(delta) < EPS) {
            if (minX[1] - EPS > inputData[i][1]) {
                minX = inputData[i];
            }
        }

        delta = maxX[0] - inputData[i][0];
        if (delta < -EPS) {
            maxX = inputData[i];
        } else if (fabs(delta) < EPS) {
            if (maxX[1] + EPS < inputData[i][1]) {
                maxX = inputData[i];
            }
        }
    }

    // R("leftEdge  " << minX[0] << ", " << minX[1]);
    // R("rightEdge " << maxX[0] << ", " << maxX[1]);
    data_t* topPlane = new data_t,
          * botPlane = new data_t;
    for (auto& pt : inputData) {
        int o = orientation(minX[0], minX[1],
                            maxX[0], maxX[1],
                            pt[0],   pt[1]);
        if (o == 1) {
            topPlane -> push_back(pt);
        } else if (o == 2) {
            botPlane -> push_back(pt);
            // R("bot " << pt[0] << ", " << pt[1]);
        } else {
            // R("col " << pt[0] << ", " << pt[1]);
        }
    }

    output.add(maxX);
    if (minX != maxX) {
        output.add(minX);
    }

    // recursive part

    recNaive(minX, maxX, topPlane);
    recNaive(maxX, minX, botPlane);

    delete topPlane;
    delete botPlane;

    return output;
}

}
