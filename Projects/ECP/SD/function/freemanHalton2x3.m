function pVal = freemanHalton2x3(obs)
% FREEMANHALTON2X3  FreemanHalton style EXACT test for a 2×3 contingency table.
%
% USAGE:
%   pVal = freemanHalton2x3(obs)
%
%   obs = [x11 x12 x13;
%          x21 x22 x23];
%
%   pVal is the two-sided (two-tailed) p-value.
%
% This function:
%   1) Ensures row sums, col sums are consistent.
%   2) Computes the probability of "obs" under the hypergeometric model.
%   3) Enumerates all feasible 2×3 tables with the same row & column sums,
%      summing probabilities of those "as or more extreme" (i.e. with probability <= prob(obs)).
%
% If row/column sums are large, enumeration can be slow. This code is for small counts.
%
% EXAMPLE:
%   obs = [ 5 2 1;
%           0 3 4];
%   pVal = freemanHalton2x3(obs)
%
% By [Your Name], [Date].

%% 1) Basic checks
if ~isequal(size(obs),[2 3])
    error('Input "obs" must be a 2x3 matrix of nonnegative integers.');
end
if any(obs(:)<0) || ~all(isfinite(obs(:))) || ~all(isreal(obs(:)))
    error('Observed counts must be real, finite, nonnegative.');
end

rowSum = sum(obs,2);   % [2 x 1]
colSum = sum(obs,1);   % [1 x 3]
N      = sum(rowSum);  % total

% check sums
if (rowSum(1)+rowSum(2))~=N || sum(colSum)~=N
    error('Inconsistent row/column sums');
end

%% 2) Precompute the log-factorial cost for speed
% Probability of a table T is:  ( rowSum(1)! rowSum(2)! colSum(1)! colSum(2)! colSum(3)! ) / N!  
%                               * 1 / ( T(1,1)! T(1,2)! T(1,3)! T(2,1)! T(2,2)! T(2,3)! )
% We'll do everything in log space using MATLAB's gammaln (where log(k!)=gammaln(k+1)).

% constant term K = sum of log(rowSum!) + sum of log(colSum!) - log(N!)
K = sum(gammaln(rowSum+1)) + sum(gammaln(colSum+1)) - gammaln(N+1);

% Observed table's log-prob
logProbObs = K ...
  - sum(gammaln(obs(1,:)+1)) ...
  - sum(gammaln(obs(2,:)+1));

pObs = exp(logProbObs);  % The observed table's probability

%% 3) Enumerate all feasible 2×3 tables with the same row sums & column sums
R1 = rowSum(1); % e.g. how many in row1
R2 = rowSum(2); % how many in row2
C1 = colSum(1);
C2 = colSum(2);
C3 = colSum(3);

% For each feasible x11, x12 in row1:
%   x13 = R1 - x11 - x12
% Then row2 is:
%   x21 = C1 - x11
%   x22 = C2 - x12
%   x23 = C3 - x13
% We must ensure all those are >=0 and sum up correctly.

pSum = 0;
pMid = 0;  % for mid p-correction if you want it
nFeasible = 0;

% define the observed table's log-prob
logProbThreshold = logProbObs;  % We'll do "two-tailed" by prob <= pObs

for x11 = 0 : C1
    if x11 > R1, break; end
    for x12 = 0 : C2
        if (x11+x12)>R1, break; end
        x13 = R1 - x11 - x12;
        if x13<0 || x13> C3, continue; end

        x21 = C1 - x11;
        x22 = C2 - x12;
        x23 = C3 - x13;

        if x21<0 || x22<0 || x23<0, continue; end
        if (x21 + x22 + x23)~=R2, continue; end

        % This table is feasible. Compute log-prob
        % table = [ x11 x12 x13;
        %           x21 x22 x23 ];
        logProb = K ...
          - ( gammaln(x11+1) + gammaln(x12+1) + gammaln(x13+1) ...
            + gammaln(x21+1) + gammaln(x22+1) + gammaln(x23+1) );

        pValTab = exp(logProb);

        nFeasible = nFeasible + 1;

        % Compare
        if logProb <= logProbThreshold
            % This table is "as or more extreme" => add its prob
            pSum = pSum + pValTab;
            if abs(logProb - logProbThreshold) < 1e-12
                % same prob => do mid-p if you want
                pMid = pMid + 0.5*pValTab; 
            end
        end
    end
end

% "Two-tailed" p-value is sum of prob of all tables with prob <= prob(obs)
% plus half the prob of ties (the mid-p approach is optional).
pValNoMid = pSum;
pValMid   = pMid + (pSum - pMid);

% Well return the classic pValNoMid or the mid pValMid. Let's do pValNoMid:
pVal = pValNoMid;

% Optionally display some info
% fprintf('FreemanHalton 2x3 enumerated %d feasible tables. pVal=%.4g\n', nFeasible, pVal);
end
