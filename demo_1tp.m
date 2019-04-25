% Demo for fitting the drift + diffusion model described in:
% Panichello MF, DePasquale B, Pillow JW, Buschman TJ.
% Error-correcting dynamics in visual working memory.
%
% Reduced model with only one set of drift and diffusion parameters (e.g
% for fitting designs with one timepoint)
%
% Dependencies:
% MATLAB R2016b or later
% Statistics and Machine Learning Toolbox

clear;

load('sampleData.mat','dat');
% fit model parameters
lowerBound = [0 0 0 0 0];
upperBound = [1 0 1 0 1];

res = dpFit(dat.report, dat.target, dat.nonTarget, ...
    dat.delayTime, dat.setSize, lowerBound, upperBound); 