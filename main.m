addpath(genpath(pwd))
data=importdata("sonar.mat");
rh=rho(data);
[Reduct_location,Location_index]=Heurstic_TMAEFS(data,rh);