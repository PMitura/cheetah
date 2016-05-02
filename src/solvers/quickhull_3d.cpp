#include "quickhull_3d.h"

namespace ch
{

Quickhull3D::Quickhull3D()
{
    name_ = "3D Quickhull";
}

Polyhedron& Quickhull3D::solve(const Points3D& input, Polyhedron& output)
{
    // edge cases
    if (input.getSize() == 0) {
        return output;
    } else if (input.getSize() <= 3) {
        Points3D plane;
        const data_t& inputData = input.getData();
        for (auto& i : inputData) {
            plane.add(i);
        }
        output.addFace(plane);
        return output;
    }

    solveSequential(input, output);
    return output;
}

Polyhedron& Quickhull3D::solveSequential(const Points3D& input,
                                         Polyhedron& output)
{
    const data_t& idata = input.getData();
    faces_.clear();

    findInitial(idata);

    return output;
}

void Quickhull3D::findInitial(const data_t& input)
{
    // find extreme points
    unsigned minpts[3], maxpts[3];
    for (int c = 0; c < 3; c++) {
        minpts[c] = 0;
        maxpts[c] = 0;
    }
    point_t pt;
    for (int c = 0; c < 3; c++) {
        for (unsigned i = 1; i < input.size(); i++) {
            pt = input[i];
            if (pt[0] > input[maxpts[c]][0]) {
                maxpts[c] = i;
            } else if (pt[0] < input[minpts[c]][0]) {
                minpts[c] = i;
            }
        }
    }

    // find first edge ab - farthest pair of axis extremes
    std::vector<unsigned> intlpts;
    double maxDist = -1;
    unsigned farc;
    for (int c = 0; c < 3; c++) {
        double d = dist3d(input[minpts[c]], input[maxpts[c]]);
        if (d > maxDist) {
            farc = c;
            maxDist = d;
        }
    }
    if (maxDist < EPS) {
        // all identical points
        return;
    }

    // find third point c - point farthest from ab
    QVertex a(input[minpts[farc]]), b(input[maxpts[farc]]);
    point_t ab = {b[0] - a[0], b[1] - a[1], b[2] - a[2]}, bc, abcNorm, abcTemp;
    double dd;
    maxDist = 0;
    int cidx = -1;
    for (unsigned i = 0; i < input.size(); i++) {
        // only relative distance
        bc = {input[i][0] - b[0], input[i][1] - b[1], input[i][2] - b[2]};
        abcTemp = perpend3d(ab, bc);
        dd = vectSqr3d(abcNorm);
        if (dd > maxDist + EPS) {
            maxDist = dd;
            cidx = i;
            abcNorm = abcTemp;
        }
    }
    if (cidx == -1) {
        // all collinear points
        return;
    }
    QVertex c(input[cidx]);

    // find fourth point d - point farthest from abc plane
    maxDist = 0;
    int didx = -1;
    double base = dot(c.crds_, abcNorm);
    for (unsigned i = 0; i < input.size(); i++) {
        dd = fabs(dot(input[i], abcNorm) - base);
        if (dd > maxDist + EPS) {
            didx = i;
            maxDist = dd;
        }
    }
    if (didx == -1) {
        // all coplanar points
        return;
    }
    QVertex d(input[didx]);

    // find out orientation of found tetrahedron and initialize half edge mesh
    // represented polyhedron accordingly
    for (int i = 0; i < 4; i++) {
        faces_.push_back(QFace());
    }
    std::vector<QVertex> used;
    if (dot(d.crds_, abcNorm) - base < -EPS) {
        used = {a, b, c};
        faces_[0].init(used);
        used = {d, b, a};
        faces_[1].init(used);
        used = {d, c, b};
        faces_[2].init(used);
        used = {d, a, c};
        faces_[3].init(used);

        faces_[1].edgeAt(1) -> pairWith(faces_[2].edgeAt(0));
        faces_[1].edgeAt(2) -> pairWith(faces_[0].edgeAt(1));
        faces_[2].edgeAt(1) -> pairWith(faces_[3].edgeAt(0));
        faces_[2].edgeAt(2) -> pairWith(faces_[0].edgeAt(2));
        faces_[3].edgeAt(1) -> pairWith(faces_[1].edgeAt(0));
        faces_[3].edgeAt(2) -> pairWith(faces_[0].edgeAt(0));
    } else {
        used = {a, c, b};
        faces_[0].init(used);
        used = {d, a, b};
        faces_[1].init(used);
        used = {d, b, c};
        faces_[2].init(used);
        used = {d, c, a};
        faces_[3].init(used);

        faces_[1].edgeAt(0) -> pairWith(faces_[2].edgeAt(1));
        faces_[1].edgeAt(2) -> pairWith(faces_[0].edgeAt(0));
        faces_[2].edgeAt(0) -> pairWith(faces_[3].edgeAt(1));
        faces_[2].edgeAt(2) -> pairWith(faces_[0].edgeAt(2));
        faces_[3].edgeAt(0) -> pairWith(faces_[1].edgeAt(1));
        faces_[3].edgeAt(2) -> pairWith(faces_[0].edgeAt(1));
    }

    // divide rest of the points under faces of tetrahedron
    // TODO
}

}
