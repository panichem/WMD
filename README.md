# Discrete attractor model of working memory dynamics

This repository contains an implementation of the discrete attractor model used in the paper "Error-correcting dynamics in visual working memory" by Panichello, DePasquale, Pillow, and Buschman (2019). 

Below are instructions for running a simple demo in which the model is fit to a sample dataset. 

1. Dependencies

    * MATLAB 2016b or later
    * Statistics and Machine Learning Toolbox
  
2. Installation guide (<1 minute)

    1. Download the repository folder to your desktop and unzip
    2. Launch MATLAB and navigate to the unzipped repository folder using the 'CurrentFolder' gui
    3. Enter the following in the MATLAB command window:
    ```matlab
                >> addpath(pwd)
    ```
3. Demo/Instructions for use
   1. To load a sample dataset and fit the discrete attractor model to this dataset, enter `>> demo` in the matlab command window. This script will take approximately 5 minutes to run. 
    
   2. The function *dpFit* in demo.m returns a structure *res*. This structure contains the negative log-likelihood of the model fit and the associated maximum likelihood parameters. To view a list of the returned parameters, enter `>> res`. The output should look something like this:
    ```
                >> res
                
                res  = 
                
                     struct with fields:
                     
                            nLL: 618.8411
                             ss: [2x1 double]
                         sigmaM: [2x1 double]
                         sigmaE: [2x1 double]
                          betaM: [2x1 double]
                          betaE: [2x1 double]
                            pgB: [2x1 double]
                            psB: [2x1 double]
                         sigmaR: 2.0745
    ```
    
    3. The value of specific parameters, (e.g., betaM for each set size) can be displayed using standard MATLAB conventions for indexing into structures, e.g.:
```
   >> res.betaM
   
   ans = 
   
        5.2467
        7.6534
```        
 4. For more information try `>> help dpFit`
