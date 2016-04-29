#include "quickhull_3d.h"

namespace ch
{

Polyhedron& Quickhull3D::solve(const Points3D& input, Polyhedron& output)
{
    // edge cases
    if (input.getSize() == 0) {
        return output;
    } else if (input.getSize() < 3) {
        Points3D plane;
        const data_t& inputData = input.getData();
        for (auto& i : inputData) {
            plane.add(i);
        }
        output.addFace(plane);
        return output;
    }

    solveNaive(input, output);
    return output;
}

Polyhedron& Quickhull3D::solveNaive(const Points3D& input, Polyhedron& output)
{
    return output;
}

}
