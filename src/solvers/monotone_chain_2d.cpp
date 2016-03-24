#include "solvers/monotone_chain_2d.h"

namespace ch
{

MonotoneChain2D::MonotoneChain2D()
{
    name_ = "Monotone Chain";
}

Points2D& MonotoneChain2D::solve(const Points2D& input, Points2D& output)
{
    return output;
}

}
