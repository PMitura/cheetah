#include "quickhull_3d.h"

namespace ch
{

Quickhull3D::Quickhull3D()
{
    EPS_LOC = 1e-6;
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

    for (auto ptr : faces_) {
        if (ptr -> state_ != QFace::OPEN) {
            continue;
        }
        Points3D facePts;
        QHalfEdge * edge = ptr -> edge_,
                  * ori  = edge;
        do {
            facePts.add(edge -> head_ -> crds_);
            edge = edge -> next_;
        } while (edge != ori);
        output.addFace(facePts);
        delete ptr;
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
    R("NEXT TOP " << assigned_.begin() -> first);
    int nxt = -1;
    double maxDist = 0, dd;
    // iterate over points assigned to that face
    for (auto& i : face -> assigned_) {
        dd = distPlanePoint(face, (*globIn_)[i]);
        if (dd > maxDist + EPS_LOC) {
            maxDist = dd;
            nxt = i;
        }
    }

    return {face, nxt};
}

void Quickhull3D::processVertex(QFace * face, unsigned index)
{
    // extract point from set of assigned
    point_t eyePoint = (*globIn_)[index];
    unassign(index);

    // find horizon of visibility from chosen point
    std::vector<QHalfEdge*> horizon;
    findHorizon(eyePoint, face, NULL, horizon);
    R("  HORIZON FOUND");

    // add faces around the horizon
    std::vector<QFace*> addedFaces;
    updateFaces(eyePoint, horizon, addedFaces);
    R("  FACES UPDATED");
 
    // first merge phase on new faces
    for (auto& face : addedFaces) {
        mergePhaseOne(face);
    }
    R("  MERGE 1 DONE");

    // second merge phase on new faces
    for (auto& face : addedFaces) {
        mergePhaseTwo(face);
    }
    R("  MERGE 2 DONE");

    // resolve newly unassigned points 
    assignOrphans(addedFaces);
    R("  ORPHANS ASSIGNED");
}

void Quickhull3D::mergePhaseOne(QFace * face)
{
    if (face -> state_ != QFace::OPEN) {
        return;
    }
    bool convex = 1;
    while (1) {
        QHalfEdge * edge = face -> edge_, * ori = edge;
        do {
            bool doMerge = 0;
            QFace * mergeWith = edge -> twin_ -> face_;
            double local = distPlanePoint(edge -> face_, 
                           edge -> twin_ -> face_ -> centroid_),
                   twin  = distPlanePoint(edge -> twin_ -> face_, 
                           edge -> face_ -> centroid_);
            if (face -> area_ > mergeWith -> area_) {
                if (local > -EPS_LOC) {
                    doMerge = 1;
                } else if (twin > -EPS_LOC) {
                    convex = 0;
                }
            } else {
                if (twin > -EPS_LOC) {
                    doMerge = 1;
                } else if (local > -EPS_LOC) {
                    convex = 0;
                }
            }
            if (doMerge) {
                R("    MERGE")
                std::vector<QFace*> merged;
                mergeFaces(face, edge, merged);
                for (auto& m : merged) {
                    reassignPoints(m, face);
                }
                continue;
            }
            edge = edge -> next_;
        } while (edge != ori);
        if (!convex) {
            face -> state_ = QFace::MERGED;
        }
        break;
    }
}

void Quickhull3D::mergePhaseTwo(QFace * face)
{
    if (face -> state_ != QFace::OPEN) {
        return;
    }
    while (1) {
        QHalfEdge * edge = face -> edge_, * ori = edge;
        do {
            bool doMerge = 0;
            double local = distPlanePoint(edge -> face_, 
                           edge -> twin_ -> face_ -> centroid_),
                   twin  = distPlanePoint(edge -> twin_ -> face_, 
                           edge -> face_ -> centroid_);
            if (local > -EPS_LOC || twin > -EPS_LOC) {
                doMerge = 1;
            }
            if (doMerge) {
                std::vector<QFace*> merged;
                mergeFaces(face, edge, merged);
                for (auto& m : merged) {
                    reassignPoints(m, face);
                }
                continue;
            }
            edge = edge -> next_;
        } while (edge != ori);
        break;
    }
}

void Quickhull3D::reassignPoints(QFace * from, QFace * to)
{
    for (auto& pt : from -> assigned_) {
        double d = distPlanePoint(to, (*globIn_)[pt]);
        if (d > EPS_LOC) {
            assign(pt, to);
        } else {
            orphans_.insert(pt);
        }
    }
    from -> assigned_.clear();
}

void Quickhull3D::mergeFaces(QFace * face, QHalfEdge * through,
        std::vector<QFace*>& merged)
{

    // close adjacent merged edge;
    QFace * toMerge = through -> twin_ -> face_;
    toMerge -> state_ = QFace::CLOSED;
    merged.push_back(toMerge);

    // update surrounding edges
    QHalfEdge * localPrev = through -> prev_,
              * twinPrev  = through -> twin_ -> prev_,
              * localNext = through -> next_,
              * twinNext  = through -> twin_ -> next_, 
              * curr;
    
    while (localPrev -> twin_ -> face_ == toMerge) {
        localPrev = localPrev -> prev_;
        twinNext  = twinNext  -> next_;
    }

    while (localNext -> twin_ -> face_ == toMerge) {
        twinPrev  = twinPrev   -> prev_;
        localNext = localNext  -> next_;
    }

    // update facial status of half-edges in merged face
    curr = twinNext;
    while (curr != twinPrev -> next_) {
        curr -> face_ = face;
        curr = curr -> next_;
    }

    // cannot use deleted edge as linked list representation
    if (through == face -> edge_) {
        face -> edge_ = localNext;
    }

    QFace * deletedA = mergeEdges(face, twinPrev, localNext),
          * deletedB = mergeEdges(face, localPrev, twinNext);
    if (deletedA) {
        merged.push_back(deletedA);
        deletedA -> state_ = QFace::MERGED;
    }
    if (deletedB) {
        merged.push_back(deletedB);
        deletedB -> state_ = QFace::MERGED;
    }

    // recalculate face properties
    face -> updateAttributes();
}

QFace * Quickhull3D::mergeEdges(QFace * face, QHalfEdge * edgeA,
        QHalfEdge * edgeB)
{
    QFace * deleted = NULL;

    if (edgeA -> twin_ -> face_ == edgeB -> twin_ -> face_) {
        // really merge into one
        QFace * twinFace = edgeA -> twin_ -> face_;
        QHalfEdge * newTwin = NULL;
        // if edgeA was head of list, move the head
        if (edgeA == face -> edge_) {
            face -> edge_ = edgeB;
        }
        if (twinFace -> getVerticesCount() == 3) {
            // delete whole face between merged edges
            deleted = twinFace;
            newTwin = edgeB -> twin_ -> prev_ -> twin_;
        } else {
            newTwin = edgeB -> twin_ -> next_;
            // change if head of list
            if (twinFace -> edge_ == newTwin -> prev_) {
                twinFace -> edge_ = newTwin;
            }
            newTwin -> prev_ = newTwin -> prev_ -> prev_;
            newTwin -> prev_ -> next_ = newTwin;
        }
        edgeB -> prev_ = edgeA -> prev_;
        edgeB -> prev_ -> next_ = edgeB;
        edgeB -> pairWith(newTwin);

        twinFace -> updateAttributes();
    } else {
        // just connect
        edgeA -> next_ = edgeB;
        edgeB -> prev_ = edgeA;
    }

    return deleted;
}

void Quickhull3D::assignOrphans(std::vector<QFace*>& candidates)
{
    // assign or trash all orphans
    for (auto& vtx : orphans_) {
        double maxDist = 0;
        QFace * chosen = NULL;

        // only check new faces
        for (auto& can : candidates) {
            if (can -> state_ != QFace::OPEN) {
                continue;
            }
            double d = distPlanePoint(can, (*globIn_)[vtx]);
            if (d > maxDist + EPS_LOC) {
                maxDist = d;
                chosen = can;
            }
        }
        if (chosen) {
            assign(vtx, chosen);
        }
        // else trashed
    }
}

void Quickhull3D::updateFaces(const point_t& eye,
        std::vector<QHalfEdge*>& horizon, std::vector<QFace*>& added)
{
    QHalfEdge * pre = NULL, * ori = NULL;
    for (auto& edge : horizon) {
        QHalfEdge * lateral = newFaceFromEdge(eye, edge);
        added.push_back(lateral -> face_);
        // check if found lateral edges of previous circular face
        if (pre == NULL) {
            ori = lateral;
        } else {
            lateral -> next_ -> pairWith(pre);
        }
        pre = lateral;
    }
    ori -> next_ -> pairWith(pre);
}

QHalfEdge * Quickhull3D::newFaceFromEdge(const point_t& eye, QHalfEdge * edge)
{
    // triangle forming the face
    std::vector<QVertex> faceVertices = {eye, *(edge -> head_),
                                         *(edge -> prev_ -> head_)};
    faces_.push_back(new QFace());
    faces_.back() -> init(faceVertices);
    QHalfEdge * faceEdge = faces_.back() -> edge_;
    faceEdge -> prev_ -> pairWith(edge -> twin_);
    return faceEdge;
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
    /*
    QHalfEdge * curr = (!through) ? on -> edge_ : through,
              * ori = curr;
    */
    QHalfEdge * curr, * ori = through;
    if (ori == NULL) {
        ori = on -> edge_;
        curr = ori;
    } else {
        curr = ori -> next_;
    }

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
    if (maxDist < EPS_LOC) {
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
        if (dd > maxDist + EPS_LOC) {
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
        if (dd > maxDist + EPS_LOC) {
            didx = i;
            maxDist = dd;
        }
    }
    if (didx == -1) {
        // all coplanar points
        return;
    }
    QVertex d(input[didx]);
    R("INITIAL " << aidx << ", " << bidx << ", " << cidx << ", " << didx);

    // find out orientation of found tetrahedron and initialize half edge mesh
    // represented polyhedron accordingly
    for (int i = 0; i < 4; i++) {
        faces_.push_back(new QFace());
    }
    std::vector<QVertex> used;
    if (dot(d.crds_, abcNorm) - base < -EPS_LOC) {
        used = {a, b, c};
        faces_[0] -> init(used);
        used = {d, b, a};
        faces_[1] -> init(used);
        used = {d, c, b};
        faces_[2] -> init(used);
        used = {d, a, c};
        faces_[3] -> init(used);

        faces_[1] -> edgeAt(1) -> pairWith(faces_[2] -> edgeAt(0));
        faces_[1] -> edgeAt(2) -> pairWith(faces_[0] -> edgeAt(1));
        faces_[2] -> edgeAt(1) -> pairWith(faces_[3] -> edgeAt(0));
        faces_[2] -> edgeAt(2) -> pairWith(faces_[0] -> edgeAt(2));
        faces_[3] -> edgeAt(1) -> pairWith(faces_[1] -> edgeAt(0));
        faces_[3] -> edgeAt(2) -> pairWith(faces_[0] -> edgeAt(0));
    } else {
        used = {a, c, b};
        faces_[0] -> init(used);
        used = {d, a, b};
        faces_[1] -> init(used);
        used = {d, b, c};
        faces_[2] -> init(used);
        used = {d, c, a};
        faces_[3] -> init(used);

        faces_[1] -> edgeAt(0) -> pairWith(faces_[2] -> edgeAt(1));
        faces_[1] -> edgeAt(2) -> pairWith(faces_[0] -> edgeAt(0));
        faces_[2] -> edgeAt(0) -> pairWith(faces_[3] -> edgeAt(1));
        faces_[2] -> edgeAt(2) -> pairWith(faces_[0] -> edgeAt(2));
        faces_[3] -> edgeAt(0) -> pairWith(faces_[1] -> edgeAt(1));
        faces_[3] -> edgeAt(2) -> pairWith(faces_[0] -> edgeAt(1));
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
            dd = distPlanePoint(faces_[f], input[i]);
            if (dd > maxDist + EPS_LOC) {
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
    assigned_[vertexID] = faces_[faceID];
    faces_[faceID] -> assigned_.insert(vertexID);
}

void Quickhull3D::assign(unsigned vertexID, QFace * face)
{
    R("ASSIGN " << vertexID);
    assigned_[vertexID] = face;
    face -> assigned_.insert(vertexID);
}

void Quickhull3D::unassign(unsigned vertexID)
{
    R("UNASSIGN " << vertexID);
    assigned_[vertexID] -> assigned_.erase(vertexID);
    assigned_.erase(vertexID);
}

}
