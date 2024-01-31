# TestOVM
Test for tetrahedral mesh subdivision with the library OpenVolumeMesh.


## Description

A tetrahedral mesh consisting of two tetrahedra sharing one face is initially generated.
Then, the application subdivides all tets incident with different vertices.


## Compilation Instructions

```shell
git submodule update --init --recursive
cmake -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build --config Debug --parallel 4
```
