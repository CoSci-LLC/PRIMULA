# PRIMULA
#### Version v1.0
*A C++ version of PRobabilistic MUltidimensional shallow Landslide Analysis*

## Dependencies
All dependencies required for this library will be download via the CMake
build process

## Building and Installation

### Linux

```bash
$ mkdir build
$ cd build 
$ cmake .. -DCMAKE_INSTALL_PREFIX=<install location>
$ make
$ make install
```

### Windows
Visual Studio 2019 is the preferred method for building PRIMULA on Windows.
CMake is also required.

Clone or download PRIMULA and open a VS 2019 Developer Console in
that directory. Then:

```bash
$ mkdir build
$ cd build
$ cmake ..
$ cmake --build .
```

## Known bugs

PRIMULA will not compile under Windows because of an OpenMP issue. To be
resolved

