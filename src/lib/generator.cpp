#include "lib/generator.h"

namespace ch
{


bool Generator2D::genUniformCircle(long long n, long long h,
                      double radius, Points2D& points)
{
    if (n < h)
        return false;

    // generate points on circle (= on hull) uniformly
    double step = (2*PI) / h;
    for (long long i = 0; i < h; i++) {
        // based on parametric equation for circle
        double px = radius * cos(i * step),
               py = radius * sin(i * step);
        points.add({px, py});
    }

    // remaining interior points
    generateInCombination(n - h, std::min(h, (long long) 50), points);

    return true;
}

void Generator2D::generateInCombination(long long n, int of, Points2D& points)
{
    srand(time(NULL));
    const data_t& hullData = points.getData();
    double comb[50];
    long long sum;
    for (long long i = 0; i < n; i++) {
        sum = 0;
        for (int j = 0; j < of; j++) {
            comb[j] = rand();
            sum += comb[j];
        }
        for (int j = 0; j < of; j++) {
            double newX = hullData[j][0] * (comb[j] / sum),
                   newY = hullData[j][1] * (comb[j] / sum);
            points.add({newX, newY});
        }
    }
}

}
