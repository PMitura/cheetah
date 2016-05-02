#pragma once

#include "lib/geometry.h"
#include "lib/structures.h"

/**
 * Structures forming Half Edge Mesh data structure
 */

namespace ch
{

struct QVertex;
struct QHalfEdge;
struct QFace;

struct QVertex {
    QVertex(point_t pt)
        : next_(NULL), prev_(NULL), edge_(NULL), crds_(pt) {}
    QVertex() : QVertex({0, 0, 0}) {}
    QVertex * next_, * prev_;
    QHalfEdge * edge_;
    point_t crds_;
    double operator[](unsigned x) {
        return crds_[x];
    }
};

struct QHalfEdge {
    QHalfEdge() : head_(NULL), twin_(NULL), prev_(NULL), next_(NULL),
        face_(NULL) {}
    QVertex * head_;
    QHalfEdge * twin_, * prev_, * next_;
    QFace * face_;

    void pairWith(QHalfEdge * twin);
};

struct QFace {
    QFace() : next_(NULL), prev_(NULL), edge_(NULL) {}
    ~QFace();
    QFace * next_, * prev_;
    QHalfEdge * edge_;
    std::vector<unsigned> conflicts_;
    point_t normal_, centroid_;
    double offset_;

    /** Initialize face as cyclic linked list of edges from vertices */
    void init(std::vector<QVertex>& vertices);

    /** Get n-th edge on a face */
    QHalfEdge * edgeAt(unsigned pos);
};

}
