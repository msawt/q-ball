# Q-Ball
  This is an implementation of the Q-Ball reconstruction algorithm (Tuch et al. 2004) in C. Specifically, it uses the "soft equator" A = ϕ(cos−1 |UTQ|). Built using tracktools' platform for calculating maxima and loading in the data.

## Usage
  Navigate to the directory that q-ball is located.
  
  ./qbi_recon -data [data] -bval [bval] -bvec [bvec] -mask [mask] -ndir [ndir] -odir [odir] -datadir [datadir] -S0 [S0]
  
  For help, use ./qbi_recon -h  
  
## Requirements
  znzlib
  
  gsl
 
 
### Written by Mason Sawtell (msawtell@uci.edu)
