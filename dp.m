function [P, xc, df] = dp(T,NT,D,SS,sigmaM,sigmaE,betaM,betaE,puB,psB,w,sigmaR) %#codegen
%DP performs density propogation in circular space

% Output:
%  P:       nBins x nTrials density matrix
%  xc:      nTrials x 1 vector of bin centers
%  df:      nBins x 1 drift function


%% Prep vars

n          = 100; 
varI       = 1e-2; 
kBase      = sd2k((2*pi/numel(w))); 
uSS        = unique(SS);
nSS        = numel(uSS);
nTrials    = numel(T);

%% Prepare density prop with matrix exponential

dx = 2*pi/n; %bin width
xe = linspace(-pi,pi,n+1); %bin edges
xc = (xe(1:n) + dx/2)'; %bin centers

cDiff = 1/(dx^2*2);  % scale factor for diffusion
cDrft = 1/(dx*2);   % scale factor for drift

% diffusion matrix
Ddff = [-2 * cDiff, cDiff, zeros(1,n-3), cDiff];
Ddff = toeplitz([Ddff(1) fliplr(Ddff(2:end))],Ddff);

%drift matrix
Dder = [0, cDrft, zeros(1,n-3), -cDrft];
Dder = toeplitz([Dder(1) fliplr(Dder(2:end))],Dder);

%drift function
muVM = -pi:(2 * pi / numel(w)):(pi - 2 * pi / numel(w));
df   = sum(w' .* (kBase * sin(muVM - xc) .* exp(kBase * cos(muVM - xc))) ...
    ./sum(2 * pi * besseli(0,kBase)),2);
df   = df / max(abs(df));

%Fokker-Planck
for iSS = 1:nSS
    M{iSS} = Ddff * sigmaM(iSS) - Dder * diag(betaM(iSS) * df);
    E{iSS} = Ddff * sigmaE(iSS) - Dder * diag(betaE(iSS) * df);
end
R = Ddff * sigmaR;

%% Mixture weights

pU  = puB(1) + SS * puB(2) + D * puB(3) + SS .* D * puB(4);
pS  = psB(1) + SS * psB(2) + D * psB(3) + SS .* D * psB(4);
pS(SS == 1) = 0;
pM = (1 - pU - pS);

%% Initialize distribution

PT = (1 / (2 * pi * besseli(0,sd2k(sqrt(varI))))) ...
    .* exp(cos(xc - T') * sd2k(sqrt(varI))) * dx;
PD = nanmean((1 / (2 * pi * besseli(0,sd2k(sqrt(varI))))) ...
    * exp(cos(xc - permute(NT,[3 1 2])) * sd2k(sqrt(varI))) * dx,3);
PD(isnan(PD)) = 0;
P0 = PT * diag(pM ./ (pM+pS)) +  PD * diag(pS ./ (pM+pS));

%% Simulate delay

%density propogation
P = nan(n,nTrials);
for iSS = 1:nSS
    tmpSS = uSS(iSS);
    for tmpD = unique(D(:)')
        ind = (D == tmpD) & (SS == tmpSS);
        P(:,ind) = expm(M{iSS} * tmpD) * expm(E{iSS}) * P0(:,ind);
    end
end
%guessing component
U = ones(n,nTrials)/n;
P = P * diag((1 - pU)) + U * diag(pU);
P = expm(R) * P;

function K = sd2k (S)
% SD2K (S)
%   Returns the Von Mises concentration parameter corresponding to
%   standard deviation S of a wrapped normal distribution.
%
%   Ref: Topics in Circular Statistics, S. R. Jammalamadaka & A. Sengupta
%
%   --> www.paulbays.com

R = exp(-S.^2/2);
K = 1./(R.^3 - 4 * R.^2 + 3 * R);
K(R < 0.85) = -0.4 + 1.39 * R(R < 0.85) + 0.43./(1 - R(R < 0.85));
K(R < 0.53) = 2 * R(R < 0.53) + R(R < 0.53).^3 + (5 * R(R < 0.53).^5)/6;

