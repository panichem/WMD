function res = dpFit(R,T,NT,D,SS,lo,hi)
%DPFIT (R,T,NT,D,SS,lo,hi)
%returns the maximum likelihood parameters of a model describing
%delayed estimation reports as the evolution of a drift and diffusion
%process in a circular space.
%
% Required Input:
% T:       nTrials x 1 vector of target values, range -pi <= X < pi
% NT:      nTrials x (max(ss)-1) matrix of non-target values, range -pi <= X < pi
% D:       nTrials x 1 vector of delay time, in seconds
% SS:      nTrials x 1 vector of set sizes
%
% Optional Input:
% lo:       5 x 1 vector of lower bounds for [sigmaM sigmaE betaM betaE sigmaR].
% hi:       5 x 1 vector of upper bounds for [sigmaM sigmaE betaM betaE sigmaR].
%
% Output:
% res.nLL:     negative log likelihood of the model
% res.ss:      unique set sizes fit
% res.sigmaM:  memory diffusion (radians/s)
% res.sigmaE:  encoding diffusion (radians/s)
% res.betaM:   memory drift (radians/s)
% res.betaE:   encoding drift (radians/s)
% res.pgB:     guessing coefficients [intercept setSize delay interaction]
% res.psB:     swap coefficients [intercept setSize delay interaction]
% res.wts:     weights for basis functions describing drift
% res.sigmaR:  decoding/response error 
%
% Ref: Panichello MF, DePasquale B, Pillow JW, Buschman TJ.
% Error-correcting dynamics in visual working memory.

rng default
nIter = 5; %number of initial search conditions

uSS = unique(SS);
nSS = numel(uSS);

if ~exist('lo','var') || isempty(lo)
    lo = zeros(1,5);
end
if ~exist('hi','var') || isempty(hi)
    hi = ones(1,5);
end

lb = [lo(1) * ones(nSS,1);  lo(2) * ones(nSS,1);  lo(3) * ones(nSS,1); ...
    lo(4) * ones(nSS,1); -.5 * ones(8,1); zeros(12,1); lo(5)];
ub = [hi(1) * ones(nSS,1);  hi(2) * ones(nSS,1);  hi(3) * ones(nSS,1); ...
    hi(4) * ones(nSS,1); .5 * ones(8,1); ones(12,1); hi(5)];
x0 = genX0(nIter,R,T,NT,D,SS,lb,ub);

maxfeval = 10000;
options = optimset('TolX',1e-8,'maxfunevals',maxfeval, ...
    'maxiter',round(maxfeval/3),'GradObj','off','Display','notify', ...
    'DerivativeCheck','off','Algorithm','interior-point');

disp('starting fitting procedure. this will take a few minutes...');
tic;

loglik = inf;
p = NaN;
for iIter = 1:nIter
    fprintf(1,'iter %d of %d\n',iIter,nIter);
    
    [x, ll] = fmincon(@(x) dpLL(R,T,NT,D,SS,x),x0(iIter,:)', ...
        [],[],[],[],lb,ub,@(x)mixCon(x,D,SS),options);
    
    if ll < loglik
        loglik = ll;
        p = x;
    end
end

fprintf(1,'elapsed time is %.1f minutes\n',toc./60)

res.nLL      = loglik; %negative log likelihood
res.ss       = uSS;
res.sigmaM   = p((0 * nSS + 1):(1 * nSS));
res.sigmaE   = p((1 * nSS + 1):(2 * nSS));
res.betaM    = p((2 * nSS + 1):(3 * nSS));
res.betaE    = p((3 * nSS + 1):(4 * nSS));
res.pgB      = p((4 * nSS + 1):(4 * nSS + 4));
res.psB      = p((4 * nSS + 5):(4 * nSS + 8));
res.wts      = p((4 * nSS + 9):end-1);
res.wts      = res.wts / sum(res.wts);
res.sigmaR   = p(end);



