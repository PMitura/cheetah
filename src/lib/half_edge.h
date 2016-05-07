#pragma once

#include <set>

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
    QVertex(point_t pt) : edge_(NULL), assignedFace_(NULL), crds_(pt) {}
    QVertex() : QVertex({0, 0, 0}) {}
    QHalfEdge * edge_;
    QFace * assignedFace_;
    point_t crds_;
    double operator[](unsigned x) {
        return crds_[x];
    }
};

struct QHalfEdge {
    QHalfEdge(QFace * face) : head_(NULL), twin_(NULL), prev_(NULL), 
                              next_(NULL), face_(face) {}
    QVertex * head_;
    QHalfEdge * twin_, * prev_, * next_;
    QFace * face_;

    void pairWith(QHalfEdge * twin);
};

struct QFace {
    QFace() : next_(NULL), prev_(NULL), edge_(NULL), state_(OPEN) {}
    ~QFace();
    QFace * next_, * prev_;
    QHalfEdge * edge_;
    std::set<unsigned> assigned_;
    point_t normal_, centroid_;
    double offset_, area_;

    enum State {OPEN, MERGED, CLOSED};
    State state_;

    /** Initialize face as cyclic linked list of edges from vertices */
    void init(std::vector<QVertex>& vertices);

    /** Computes offset, area and cetroid */
    void updateAttributes();

    /** Returns number of vertices on perimeter */
    int getVerticesCount();

    /** Get n-th edge on a face */
    QHalfEdge * edgeAt(unsigned pos);
};

// Half-Edge Mesh related geometry methods ------------------------------------

inline double distPlanePoint(const QFace * plane, const point_t& pt)
{
    return   plane -> normal_[0] * pt[0]
           + plane -> normal_[1] * pt[1]
           + plane -> normal_[2] * pt[2]
           - plane -> offset_;
}

inline bool visible(const QFace * plane, const point_t& pt)
{
    return distPlanePoint(plane, pt) > 1e-6;
}

}
