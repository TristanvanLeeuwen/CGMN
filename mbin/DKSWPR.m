function x = DKSWPR(R,idx,x,b,w)
% double Kaczmarz sweep for band-storage matrices. Calls sweepR_mex.c
%
% use:
%   x = DKSWPR(R,idx,x,b,w)
%
% input:
%   R   - matrix bands, see mat2R.
%   idx - column indices, see mat2R
%   x   - input vector 
%   b   - r.h.s. vector 
%   w   - relaxation parameter
%
% Tristan van Leeuwen, 2013
% tristan.vanleeuwen@gmail.com

N = size(R,1);

if isempty(x)
    x = zeros(N,1);
end
if isempty(b)
    b = zeros(N,1);
end

% forward sweep
x = sweepR_mex(R,idx,x,b,w,1);

% backward sweep
x = sweepR_mex(R,idx,x,b,w,-1);
