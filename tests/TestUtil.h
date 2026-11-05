#pragma once
#include <cstdio>
#include <cstdlib>

static int g_failures = 0;
#define CHECK(cond)                                                                        \
    do {                                                                                   \
        if (!(cond)) { std::printf("  FAILED %s:%d  %s\n", __FILE__, __LINE__, #cond); ++g_failures; } \
    } while (0)
#define TEST_END() do { std::printf("%s\n", g_failures ? "TEST FAILED" : "TEST PASSED"); return g_failures ? 1 : 0; } while (0)
