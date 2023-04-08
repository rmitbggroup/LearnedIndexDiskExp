Build Runnable Program
- Build ALEX
./build.sh => build/benchmark

- Build b+-tree and fitting
g++ -O3 test2.cpp -std=c++17 -I./stx-btree-0.9/include -o test

- Build LIPP
g++ -O3 -std=c++17 test.cpp -o test

-Build PGM
g++ -O3 examples/updates.cpp -std=c++17 -I./include -o test

DATASET
- SOSD
https://github.com/learnedsystems/SOSD/blob/master/scripts/download.sh
- GRE
https://github.com/gre4index/GRE/blob/master/datasets/download.sh

