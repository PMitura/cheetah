#pragma once

namespace ch
{

/**
 * Wrapper for running performance tests of implemented algorithms
 */
class PerfTest
{
    public:
        /** runs all scheduled tests */
        void runAllTests();

        bool runGeneratedTest(int n, int h, int dim);

    private:

};

}
