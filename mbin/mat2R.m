function [R,idx] = mat2R(A)
% convert sparse matrix to band storate format
%
% use:
%   [R,idx] = mat2R(A);
%   
% input:
%   A - sparse matrix
%   
% output:
%   R   - bands of the matrix
%   idx - indices such that R(i,j) = A(i,i+idx(j))
%
% Tristan van Leeuwen, 2013
% tristan.vanleeuwen@gmail.com


%% get just the rows:
[R,idx] = spdiags(A);
for k = 1:length(idx)
    m = size(A,2);
    i = idx(k);
    R(:,k) = [zeros(-i,1);R(max(1+i,1):min(m+i,m),k);zeros(i,1)];
end