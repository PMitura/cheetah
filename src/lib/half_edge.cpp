#include "lib/half_edge.h"

namespace ch
{

void QHalfEdge::pairWith(QHalfEdge * twin)
{
    twin_ = twin;
    twin -> twin_ = this;
}

QFace::~QFace()
{
    QHalfEdge * curr, * prev;
    if (edge_ != NULL) {
        curr = edge_;
        curr -> prev_ -> next_ = NULL;
        while (curr -> next_ != NULL) {
            delete curr -> head_;
            prev = curr;
            curr = curr -> next_;
            delete prev;
        }
        delete curr -> head_;
        delete curr;
    }
}

void QFace::init(std::vector<QVertex>& vertices)
{
    // cannot initalize with less than 3 vertices
    if (vertices.size() < 3) {
        R("Failed face initialization")
        return; // shouldn't be ever called
    }

    // initialize edges in a circular linked list
    QHalfEdge * prev = NULL;
    for (auto& v : vertices) {
        QHalfEdge * curr = new QHalfEdge();
        curr -> face_ = this;
        curr -> head_ = new QVertex(v.crds_);
        curr -> head_ -> edge_ = curr;
        if (prev == NULL) {
            edge_ = curr;
        } else {
            prev -> next_ = curr;
            curr -> prev_ = prev;
        }
        prev = curr;
    }
    prev -> next_ = edge_;
    edge_ -> prev_ = prev;

    // find normal vector of face plane
    QHalfEdge * e1 = edge_ -> next_, * e2 = e1 -> next_;
    point_t pt0 = edge_ -> head_ -> crds_,
            pt2 = e1    -> head_ -> crds_,
            v02 = {pt2[0] - pt0[0], pt2[1] - pt0[1], pt2[2] - pt0[2]},
            v01, prp;
    normal_ = {0, 0, 0};
    while (e2 != edge_) {
        v01 = v02;
        pt2 = e2 -> head_ -> crds_;
        v02 = {pt2[0] - pt0[0], pt2[1] - pt0[1], pt2[2] - pt0[2]};
        prp = perpend3d(v01, v02);
        normal_[0] += prp[0];
        normal_[1] += prp[1];
        normal_[2] += prp[2];

        e1 = e2;
        e2 = e2 -> next_;
    }
    double scaleFactor = vectLen3d(normal_);
    normal_[0] /= scaleFactor;
    normal_[1] /= scaleFactor;
    normal_[2] /= scaleFactor;

    // find centroid of face
    centroid_ = {0, 0, 0};
    QHalfEdge * curr = edge_;
    do {
        centroid_[0] += curr -> head_ -> crds_[0];
        centroid_[1] += curr -> head_ -> crds_[1];
        centroid_[2] += curr -> head_ -> crds_[2];
        curr = curr -> next_;
    } while (curr != edge_);
    centroid_[0] /= vertices.size();
    centroid_[1] /= vertices.size();
    centroid_[2] /= vertices.size();

    // set plane dot product offset
    offset_ = dot(normal_, centroid_);
}

QHalfEdge * QFace::edgeAt(unsigned pos)
{
    QHalfEdge * e = edge_;
    for (unsigned i = 0; i < pos; i++) {
        e = e -> next_;
    }
    return e;
}

}
