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
    globIn_ = &idata;
    faces_.clear();
    assigned_.clear();
    orphans_.clear();

    findInitial(idata);
    // iterate over faces until done
    while (1) {
        std::pair<QFace*, int> eyePoint = nextVertex();
        if (eyePoint.second == -1) {
            break;
        }
        processVertex(eyePoint.first, eyePoint.second);
    }

    return output;
}

std::pair<QFace*, int> Quickhull3D::nextVertex()
{
    if (assigned_.empty()) {
        return {NULL, -1};
    }
    // find face of top point
    QFace * face = assigned_.begin() -> second;
    int nxt = -1;
    double maxDist = 0, dd;
    // iterate over points assigned to that face
    for (auto& i : face -> assigned_) {
        dd = distPlanePoint(face, (*globIn_)[i]);
        if (dd > maxDist + EPS) {
            maxDist = dd;
            nxt = i;
        }
    }

    return {face, nxt};
}

void Quickhull3D::processVertex(QFace * face, unsigned index)
{
    point_t eyePoint = (*globIn_)[index];
    unassign(index);

    std::vector<QHalfEdge*> horizon;
    findHorizon(eyePoint, face, NULL, horizon);

    std::vector<QFace*> addedFaces;
    // TODO updateFaces(eyePoint, horizon, addedFaces);
}

void Quickhull3D::findHorizon(const point_t& of, QFace * on, QHalfEdge * through,
                              std::vector<QHalfEdge*>& horizon)
{
    // free points from face we are looking at and close it
    for (auto& i : on -> assigned_) {
        orphans_.insert(i);
    }
    on -> assigned_.clear();
    on -> state_ = QFace::CLOSED;

    // pick next edge
    QHalfEdge * curr = (!through) ? on -> edge_ : through,
              * ori = curr;

    // DFS iterate over edges
    do {
        QHalfEdge * twin = curr -> twin_;
        QFace * twinFace = twin -> face_;
        if (twinFace -> state_ == QFace::OPEN) {
            if (visible(twinFace, of)) {
                findHorizon(of, twinFace, twin, horizon);
            } else {
                horizon.push_back(curr);
            }
        }
        curr = curr -> next_;
    } while (curr != ori);
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
    int aidx = minpts[farc], bidx = maxpts[farc];
    QVertex a(input[aidx]), b(input[bidx]);
    point_t ab = {b[0] - a[0], b[1] - a[1], b[2] - a[2]}, bc, abcNorm, abcTemp;

    // find third point c - point farthest from ab
    double dd;
    maxDist = 0;
    int cidx = -1;
    for (unsigned i = 0; i < input.size(); i++) {
        // only relative distance
        bc = {input[i][0] - b[0], input[i][1] - b[1], input[i][2] - b[2]};
        abcTemp = perpendNormal3d(ab, bc);
        dd = vectSqr3d(abcTemp);
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
    point_t curr;
    for (int i = 0; i < (int) input.size(); i++) {
        // do not reuse tetrahedron vertices
        if (i == aidx || i == bidx || i == cidx || i == didx) {
            continue;
        }

        maxDist = 0;
        int faceID = -1;
        for (int f = 0; f < 4; f++) {
            dd = distPlanePoint(&(faces_[f]), input[i]);
            if (dd > maxDist + EPS) {
                faceID = f;
                maxDist = dd;
            }
        }
        if (faceID >= 0) {
            assign(i, faceID);
        }
    }
}

void Quickhull3D::assign(unsigned vertexID, unsigned faceID)
{
    R("ASSIGN " << vertexID);
    assigned_[vertexID] = &(faces_[faceID]);
    faces_[faceID].assigned_.insert(vertexID);
}

void Quickhull3D::unassign(unsigned vertexID)
{
    R("UNASSIGN " << vertexID);
    assigned_[vertexID] -> assigned_.erase(vertexID);
    assigned_.erase(vertexID);
}

}
