#pragma once

#include <omp.h>

#include "approximators/bfp2d.h"

#include "lib/structures.h"

#include "solvers/chan_2d.h"
#include "solvers/graham_scan_2d.h"
#include "solvers/jarvis_scan_2d.h"
#include "solvers/jarvis_scan_3d.h"
#include "solvers/monotone_chain_2d.h"
#include "solvers/quickhull_2d.h"

/**
 * Header only interfaces for solvers, used in library form of ditribution.
 */

namespace ch
{

/**
 * Uses Quickhull algorithm to find convex hull of given 2D point set
 *
 * @param input Input set of points
 * @param output Reference to set of points containing convex hull
 *
 * @return Same as output param, reference to result
 */
Points2D& findHull(const Points2D& input, Points2D& output);

/** Enumerates types of 2d solvers */
enum SolverType {JARVIS, GRAHAM, ANDREW, QUICKHULL, CHAN};

/**
 * Uses chosen algorithm to find convex hull of given point set
 *
 * @param input Input set of points
 * @param output Reference to set of points containing convex hull
 * @param type SolverType of chosen algorithm
 *
 * @return Same as output param, reference to result
 */
Points2D& findHull(const Points2D& input, Points2D& output, SolverType type);

/** Parallel version of findHull, allows to choose number of threads */
Points2D& findHullParallel(const Points2D& input, Points2D& output, int thr);

/**
 * Parallel version of findHull with algorithm selection,
 * allows to choose number of threads 
 */
Points2D& findHullParallel(const Points2D& input, Points2D& output, 
        SolverType type, int thr);


/**
 * Approximates convex hull of given set of points using BFP approximation
 *
 * @param input Input set of points
 * @param output Approximation of convex hull
 *
 * @return Same as output
 */
Points2D& approximateHull(const Points2D& input, Points2D& output);

/**
 * Uses Jarvis algorithm to find convex hull of given set of 3D points
 *
 * @param input Input set of points
 * @param output Convex hull of input, represented by list of facets
 *
 * @return same as output
 */
Polyhedron& findHull3D(const Points3D input, Polyhedron& output);

}
