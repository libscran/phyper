# Hypergeometric tail calculations

![Unit tests](https://github.com/libscran/phyper/actions/workflows/run-tests.yaml/badge.svg)
![Documentation](https://github.com/libscran/phyper/actions/workflows/doxygenate.yaml/badge.svg)
[![Codecov](https://codecov.io/gh/libscran/phyper/graph/badge.svg?token=5vjoQw0cjH)](https://codecov.io/gh/libscran/phyper)

## Overview

The **phyper** library (pronounced "viper") computes the cumulative probability for a hypergeometric distribution,
intended for testing overrepresentation of pathways/signatures within the set of marker genes.
The interface is based on R's `stats::phyper()` function, with support for (log-transformed) probabilities of either tail.

## Quick start

```cpp
#include "phyper/phyper.hpp"

phyper::Options opt;
opt.upper_tail = true;
opt.log = false;

// Equivalent to stats::phyper(4, 20, 10000, 100, lower.tail=FALSE)
phyper::compute(
    /* number of marker genes from the pathway */ 5, 
    /* number of genes in the pathway */ 20, 
    /* number of genes not in the pathway */ 10000,
    /* number of marker genes */ 100,
    opt
);
```

Note that the upper-tailed cumulative probability returned by `phyper::compute()` includes the probability mass of the observed number of marker genes in the pathway.
This means that it can be directly used as the overrepresentation p-value for the pathway.
For comparable results from `stats::phyper()`, users should subtract 1 from the `q=` argument..

Check out the [reference documentation](https://libscran.github.io/phyper) for more details.

## Building projects

### CMake with `FetchContent`

If you're using CMake, you just need to add something like this to your `CMakeLists.txt`:

```cmake
include(FetchContent)

FetchContent_Declare(
  phyper
  GIT_REPOSITORY https://github.com/libscran/phyper
  GIT_TAG master # or any version of interest
)

FetchContent_MakeAvailable(phyper)
```

Then you can link to **phyper** to make the headers available during compilation:

```cmake
# For executables:
target_link_libraries(myexe libscran::phyper)

# For libaries
target_link_libraries(mylib INTERFACE libscran::phyper)
```

### CMake with `find_package()`

```cmake
find_package(libscran_phyper CONFIG REQUIRED)
target_link_libraries(mylib INTERFACE libscran::phyper)
```

To install the library, use:

```sh
mkdir build && cd build
cmake .. -DNENESUB_TESTS=OFF
cmake --build . --target install
```

### Manual

If you're not using CMake, the simple approach is to just copy the files in `include/` - either directly or with Git submodules - and include their path during compilation with, e.g., GCC's `-I`.T
