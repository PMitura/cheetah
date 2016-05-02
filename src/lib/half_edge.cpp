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
        while (curr -> next_ != NULL) {
            edge_ -> prev_ -> next_ = NULL;
            prev = curr;
            curr = curr -> next_;
            delete prev -> tail_;
            delete prev;
        }
    }
}

void QFace::init(std::vector<QVertex>& vertices)
{
    QHalfEdge * prev = NULL;
    for (auto& v : vertices) {
        QHalfEdge * curr = new QHalfEdge();
        curr -> face_ = this;
        curr -> tail_ = new QVertex(v);
        if (prev == NULL) {
            edge_ = curr;
        } else {
            prev -> next_ = curr;
            curr -> prev_ = prev;
        }
        prev = curr;
    }
    if (edge_ != NULL) {
        prev -> next_ = edge_;
        edge_ -> prev_ = prev;
    }
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
