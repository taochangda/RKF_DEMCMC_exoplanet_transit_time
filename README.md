Version 1.0
Tao Changda, September, 2019, All rights reserved.
EMAIL:risheng828@126.com


This is the programs trying to  caculate to transit time of exoplanets , based on DE-MCMC and RKF algorithms

1、Notes for program or script:
1),"import_source_data.m" : import observed data, set range of fitting parameter
2), "for_initialization_RKF78.m" or "for_initialization_RKF45.m":  initialization and pre-process
3), "main_7Np_RKF78_20190301.m" or "main_7Np_RKF45_20190301.m": main program for  fitting parameters using Differential Evolution (DE) algorithm,DE-MCMC(Markov Chain Monte Carlo) algorithm
4), "Nbody_C_7Np_RKF78.c" or "Nbody_C_7Np_RKF45.c": soure C code for searching transit time using RKF78 or RKF45. English and Chinese are noted in "Nbody_C_7Np_RKF78.c",but Chinese only or is not noted in "Nbody_C_7Np_RKF45.c".

2、Steps for build running environment and running scripts 
1),Build Parallel Computing Cluster for parallel for loop（parfor),refer to matlab offical guide for more help

2)，Build MEX function from C code(refer to matlab offical guide for more help about "mex")

clear Nbody_C_7Np_RKF45.mexw64
delete Nbody_C_7Np_RKF45.mexw64
mex Nbody_C_7Np_RKF45.c

clear Nbody_C_7Np_RKF78.mexw64
delete Nbody_C_7Np_RKF78.mexw64
mex Nbody_C_7Np_RKF78.c

3),copy Nbody_C_7Np_RKF78.mexw64 and Nbody_C_7Np_RKF45.mexw64 to all the computing Nodes(including client and workers) in the correct paths
Notes: pls clear and delete the old flies(if there were) successfully in Nodes(including client and workers) before copy new files,

4),set just the client machine for runing "main_7Np_RKF78_20190301.m" or "main_7Np_RKF45_20190301.m" automatically after reboot but not log in. and "main_startup.bat" is an example(there is no ".m" in the runnig file name)

5),run the script

firstly, run "import_source_data.m"
secondly, run "for_initialization_RKF78.m" or "for_initialization_RKF45.m"
thirdly, run " main_7Np_RKF78_20190301.m" or "main_7Np_RKF45_20190301.m"

lastly run "result_process.m" if correct result is get
