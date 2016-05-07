#pragma once

#include <vector>
#include <map>

#include "lib/geometry.h"
#include "lib/half_edge.h"
#include "lib/structures.h"
#include "solvers/solver_3d.h"

namespace ch
{

/**
 * Quickhull modified for convex hulls in 3D
 */
class Quickhull3D : public Solver3D
{
    public:
        Quickhull3D();
        /**
         * Find convex hull of input data set.
         *
         * @param inputSet Input data set
         * @return Convex hull of input data set.
         */
        Polyhedron& solve(const Points3D& input, Polyhedron& output);

    private:
        /** Beginning point for sequential solving */
        Polyhedron& solveSequential(const Points3D& input, Polyhedron& output);

        /** Gives next vertex to be processed */
        std::pair<QFace*, int> nextVertex();

        /** Processes given vertex, adding it to partial hull */
        void processVertex(QFace * face, unsigned index);

        /** Merge newly made face if needed */
        void mergePhaseOne(QFace * face);

        /** Merge newly made face if needed, second pass */
        void mergePhaseTwo(QFace * face);

        /** Reassign freed points to face or orphan them */
        void reassignPoints(QFace * from, QFace * to);

        /** Merge face through an edge */
        void mergeFaces(QFace * face, QHalfEdge * through,
                std::vector<QFace*>& merged);

        /** Merge adjacent edges on a face, return deleted face if needed */
        QFace * mergeEdges(QFace * face, QHalfEdge * edgeA, QHalfEdge * edgeB);

        /** Assign orphaned points to new candidate faces */
        void assignOrphans(std::vector<QFace*>& candidates);

        /** Updates faces according to found horizon */
        void updateFaces(const point_t& eye, std::vector<QHalfEdge*>& horizon,
                std::vector<QFace*>& added);

        /** creates new face from given eyepoint and horizon edge */
        QHalfEdge * newFaceFromEdge(const point_t& eye, QHalfEdge * edge);

        /**
         * Finds closed list of edges forming horizon visible from given point
         */
        void findHorizon(const point_t& of, QFace * on, QHalfEdge * through,
                         std::vector<QHalfEdge*>& horizon);
        /** 
         * Finds initial tetrahedron.
         *
         * @param input all input pts
         * @return tetrahedron, or empty polyhedron in case of degenerate input
         */
        void findInitial(const data_t& input);

        /** Assigns given vertice to given face (overrides for id and ptr) */
        void assign(unsigned vertexID, unsigned faceID);
        void assign(unsigned vertexID, QFace * face);

        /** Opposite of assign */
        void unassign(unsigned vertexID);

        /** Input data set made global */
        const data_t* globIn_;

        /** Global list of found faces in half-edge mesh structure */
        std::vector<QFace*> faces_;

        /** Global list of points assigned to some face */
        std::map<unsigned, QFace*> assigned_;

        /** Unassigned points */
        std::set<unsigned> orphans_;

        /** Local EPS value */
        double EPS_LOC;
};

}
