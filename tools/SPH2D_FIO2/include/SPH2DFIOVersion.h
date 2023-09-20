#pragma once

#define SPH2D_FIO_VERSION_MAJOR 2
#define SPH2D_FIO_VERSION_MINOR 3
#define SPH2D_FIO_VERSION_PATCH 0

// 2.2.0 - start version control
// 2.2.1 - fix loading layers with maxn greater than default
// 2.3.0 - fix opening experiment with no data directory
//       - filesystem::path in interface instead of string
//       - if data dir is empty then try to start from dump dir
//       - TimeLayers are constructed from ntotal instead of maxn
