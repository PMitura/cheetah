#include "lib/generator.h"

namespace ch
{


bool Generator2D::genUniformCircle(long long n, long long h,
                      double radius, Points2D& points)
{
    if (n < h || n <= 0 || radius < EPS || h <= 1) {
        return false;
    }

    // generate points on circle (= on hull) uniformly
    double step = (2*PI) / h;
    for (long long i = 0; i < h; i++) {
        // based on parametric equation for circle
        double px = radius * cos(i * step),
               py = radius * sin(i * step);
        points.add({px, py});
    }

    // remaining interior points
    // generateInCombination(n - h, std::min(h, (long long) 50), points);
    generateTriangulated(n - h, std::min(h, (long long) 50), points);

    return true;
}

bool Generator2D::genRandomCircle(long long n, long long h,
                                  double radius, Points2D& points)
{
    if (n < h || n <= 0 || radius < EPS || h <= 1) {
        return false;
    }

    // generate points on circle (= on hull) randomly
    srand(time(NULL));
    std::vector<double> coefs;
    for (long long i = 0; i < h; i++) {
        coefs.push_back(((double) rand() / RAND_MAX) * (2*PI));
    }
    std::sort(coefs.begin(), coefs.end());
    for (long long i = 0; i < h; i++) {
        double px = radius * cos(coefs[i]),
               py = radius * sin(coefs[i]);
        points.add({px, py});
    }

    // remaining interior points
    generateTriangulated(n - h, std::min(h, (long long) 50), points);

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
        double newX = 0, newY = 0;
        for (int j = 0; j < of; j++) {
            newX += hullData[j][0] * (comb[j] / sum);
            newY += hullData[j][1] * (comb[j] / sum);
        }
        points.add({newX, newY});
    }
}

void Generator2D::generateTriangulated(long long n, int of, Points2D& points)
{
    if (of <= 2) {
        generateInCombination(n, of, points);
        return;
    }
    // use first point a) a pivot, diagonals are lines from pivot to all
    // remaining points
    const data_t& hullData = points.getData();

    // compute weights of triangles
    std::vector<std::pair<double, int>> weights;
    double totalArea = 0, area, cumul = 0;
    for (int i = 1; i < of - 1; i++) {
        area = fabs((hullData[0  ][0] * (hullData[i  ][1] - hullData[i+1][1]) +
                     hullData[i  ][0] * (hullData[i+1][1] - hullData[0  ][1]) +
                     hullData[i+1][0] * (hullData[0  ][1] - hullData[i  ][1])
                    ) / 2);
        totalArea += area;
        weights.push_back({area, i});
    }
    for (auto& i : weights) {
        i.first /= totalArea;
        cumul += i.first;
        i.first = cumul;
    }

    // generate points
    srand(time(NULL));
    for (int x = 0; x < n; x++) {

        // choose random triangle
        double rng = ((double) rand() / RAND_MAX);
        std::pair<double, int> cmp = {rng, 0};
        int rngIdx = std::lower_bound(weights.begin(), weights.end(), cmp)
                     -> second;

        // generate random point inside that triangle (ABC)
        double ax = hullData[0         ][0], ay = hullData[0         ][1],
               bx = hullData[rngIdx    ][0], by = hullData[rngIdx    ][1],
               cx = hullData[rngIdx + 1][0], cy = hullData[rngIdx + 1][1];

        double rng1 = ((double) rand() / RAND_MAX),
               rng2 = ((double) rand() / RAND_MAX);

        double newPtX =   (1 - sqrt(rng1)) * ax
                        + sqrt(rng1) * (1 - rng2) * bx
                        + sqrt(rng1) * rng2 * cx,
               newPtY =   (1 - sqrt(rng1)) * ay
                        + sqrt(rng1) * (1 - rng2) * by
                        + sqrt(rng1) * rng2 * cy;

        points.add({newPtX, newPtY});
    }
}

}
