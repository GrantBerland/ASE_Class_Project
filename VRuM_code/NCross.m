function [NV, LenV, V] = NCross(L1, L2)
% Normalized vector orthogonal to L1 and L2
% [NV, LenV, V] = NCross(L1, L2)
% INPUT:
%   L1, L2: [1 x 3] or [T x 3] double arrays. The size of the 1st
%           dimension of L1 and L2 must match or one must have 1 row.
%           L1 and L2 needn't be normalized.
% OUTPUT:
%   NV:     [T x 3] array, normalized cross product between rows of the
%           input matrices:
%   LenV:   [T x 1] array, length of the cross products.
%   V:      [T x 3] array, cross vector before normalization.
%
% MATHEMATICAL DEFINITION:
%   V    = cross(L1, L2)
%   LenV = norm(V, 2)
%   NV   = V / LenV
%
% If input vectors are almost parallel, the normalization will be
% weak, therefore all vectors with length smaller than SQRT(EPS) are set to Inf.
%
% Tested: Matlab 6.5, 7.7, 7.8, 7.13, WinXP/32, Win7/64
% Author: Jan Simon, Heidelberg, (C) 2009-2013  
% $License: BSD $

if ~(isequal(2, ndims(L1), ndims(L2)) && ...
     (isequal(size(L1), size(L2)) || ...
     size(L1, 1) == 1 || size(L2, 1) == 1))
   error('JSimon:NCross:BadInputSize', ...
         'NCross needs matching dimensions for L1 and L2!');
end

%smallVal = 100.0 * eps;
smallVal = 1.4901161193847656e-008;  % SQRT(EPS)
infVal   = Inf;

V = [ ( L1(:, 2) .* L2(:, 3) - L1(:, 3) .* L2(:, 2) ), ...
      ( L1(:, 3) .* L2(:, 1) - L1(:, 1) .* L2(:, 3) ), ...
      ( L1(:, 1) .* L2(:, 2) - L1(:, 2) .* L2(:, 1) ) ];

LenV = sqrt(V(:, 1).*V(:, 1) + V(:, 2).*V(:, 2) + V(:, 3).*V(:, 3));

lowI       = (LenV < smallVal);
LenV(lowI) = 1;                        % Prevent division by zero
NV         = V ./ LenV(:, [1, 1, 1]);  % Or BSXFUN of course
LenV(lowI) = 0;
badI       = (isfinite(LenV) == 0);
LenV(badI) = infVal;

% Set weak vectors to [Inf, Inf, Inf]:
weakI        = (lowI | badI);
NV(weakI, :) = infVal;
V(weakI,  :) = infVal;