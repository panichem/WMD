% Demo for fitting the drift + diffusion model described in:
% Panichello MF, DePasquale B, Pillow JW, Buschman TJ.
% Error-correcting dynamics in visual working memory.
%
% Dependencies:
% MATLAB R2016b or later
% Statistics and Machine Learning Toolbox

clear;

load('sampleData.mat','dat');
% fit model parameters
res = dpFit(dat.report, dat.target, dat.nonTarget, ...
    dat.delayTime, dat.setSize); 