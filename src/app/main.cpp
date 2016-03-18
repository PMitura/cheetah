#include <iostream>

#include "app/perftest.h"
#include "lib/structures.h"

using namespace std;

#ifndef TEST_MAIN

/** Temporary purpose as performance testbench */
int main()
{
    ch::PerfTest perfTest;
    perfTest.runAllTests();
    return 0;
}

#endif
