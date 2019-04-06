function ll = dpLL(R,T,NT,D,SS,x)
           
nSS = numel(unique(SS));
sigmaM = x((0 * nSS + 1):(1 * nSS));
sigmaE = x((1 * nSS + 1):(2 * nSS));
betaM  = x((2 * nSS + 1):(3 * nSS));
betaE  = x((3 * nSS + 1):(4 * nSS));
puB    = x((4 * nSS + 1):(4 * nSS + 4));
psB    = x((4 * nSS + 5):(4 * nSS + 8));
w      = x((4 * nSS + 9):end-1);
sigmaR = x(end);

[P, xc] = dp(T, NT, D, SS, sigmaM, sigmaE, betaM, betaE, puB, psB, w, sigmaR);
[~, i] = min(abs(xc - R'),[],1); 
ind = sub2ind(size(P),i,1:numel(T)); 
p = P(ind);
ll = sum(-log(p));
if ~isreal(ll); ll = inf; end

