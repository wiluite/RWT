# RWT
The acronym RWT stands for Ricker Wavelet Transwarp, which means that this is an example of how to do competitive wavelet transforms through task-based parallelism provided by the Transwarp package (by Christian Blume).

transform_uni: 1516
transform_omp: 8124
transform_tw:  203

Current state requires some Boost, GCC compiler 7.2+ with extensions (to compare with naive openmp version introduced as well). 
As-is. "Windows only" for now (looking forward). Without any sanity checks of the task-based version so far.


