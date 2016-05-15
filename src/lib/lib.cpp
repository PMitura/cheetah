#include "lib/lib.h"

namespace ch
{

Points2D& findHull(const Points2D& input, Points2D& output)
{
    Quickhull2D quick;
    return quick.solve(input, output);
}

Points2D& findHull(const Points2D& input, Points2D& output, SolverType type)
{
    Solver2D * solver = NULL;
    switch (type) {
        case JARVIS:
            solver = new JarvisScan2D();
            break;
        case GRAHAM:
            solver = new GrahamScan2D();
            break;
        case QUICKHULL:
            solver = new Quickhull2D();
            break;
        case CHAN:
            solver = new Chan2D();
            break;
        case ANDREW:
            solver = new MonotoneChain2D();
            break;
        default:
            solver = new Quickhull2D();
            break;
    }
    Points2D& result = solver -> solve(input, output);
    delete solver;
    return result;
}

Points2D& findHullParallel(const Points2D& input, Points2D& output, int thr)
{
    Quickhull2D quick(Quickhull2D::PARA);
    if (thr > MAX_NUM_THREADS) {
        thr = MAX_NUM_THREADS;
    }
    omp_set_num_threads(thr);
    return quick.solve(input, output);
}

Points2D& findHullParallel(const Points2D& input, Points2D& output,
        SolverType type, int thr)
{
    Solver2D * solver = NULL;
    switch (type) {
        case JARVIS:
            solver = new JarvisScan2D(JarvisScan2D::PARA_DOUBLE);
            break;
        case GRAHAM:
            solver = new GrahamScan2D(GrahamScan2D::PARA);
            break;
        case QUICKHULL:
            solver = new Quickhull2D(Quickhull2D::PARA);
            break;
        case CHAN:
            solver = new Chan2D(Chan2D::PARA_COMBO);
            break;
        case ANDREW:
            solver = new MonotoneChain2D();
            break;
        default:
            solver = new Quickhull2D(Quickhull2D::PARA);
            break;
    }
    if (thr > MAX_NUM_THREADS) {
        thr = MAX_NUM_THREADS;
    }
    omp_set_num_threads(thr);
    Points2D& result = solver -> solve(input, output);
    delete solver;
    return result;
}

Points2D& approximateHull(const Points2D& input, Points2D& output)
{
    BFP2D bfp;
    return bfp.approximate(input, output);
}

Polyhedron& findHull3D(const Points3D input, Polyhedron& output)
{
    JarvisScan3D solver;
    return solver.solve(input, output);
}

}
