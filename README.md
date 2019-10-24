# Experiemental-design
Maintainer: Jia Liu (jia.liu@fmi.fi)

# Reference
If you use Experimental-design, please use the following reference:

Jia Liu, Jarno Vanhatalo. Bayesian model based spatiotemporal survey designs and partially observed log Gaussian Cox process. 
Journal of Spatial data analysis. http://authors.elsevier.com/authorforms/SPASTA100392/73ebeb5a162320e90dac4f3cc871c3b6

# Introduction

Experimental-design is a code toolbox in the MATLAB enviroment to do the spatial and spatiotemporal design sampling with and without the rejection algorithm.
The rejection design proposed in the above paper can give more weights
on time and/or space  which are a priori expected to be more informative through an 
intensity function of the point process.    


# Installing the toolbox

Clone (or download the source doce)  this “Experimental-design” repository and add the folder(s) to your Matlab path
or run install_design.m in the root directory of Experimental-design.



# Dependence toolbox
If you want to evaluate the design under the Gaussian process models, a dependency toolbox need to be installed, 
please see https://research.cs.aalto.fi/pml/software/gpstuff/  for more details inclduing right contact person(s), License
and etc.


# A short user guide
1. To sample the design, simply run demo_design.m in the 'design' folder. 

2. A loss/utility function is intruduced for evaluating the design with the GP models, example code files are in the 'gp_design' folder.
There is a demo,  demo_evaluate_design.m,  in the folder.

3. Evaluating design under the Gaussian process model required a dependence toolbox, gpstuff.
  

   

