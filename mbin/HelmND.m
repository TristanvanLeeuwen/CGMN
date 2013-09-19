function H = HelmND(k,d,n,bc)
% set up Helmholtz matrix in 1D, 2D or 3D with Dirichlet or absorbing bc's
%
% usage
%   H = HelmND(k,o,d,n,nb)
%
% input:
%   k       - gridded wavenumber k =  2*pi*f/c [m]
%   {o,d,n} - grid: x1 = d(1)*[0:n(1)-1]; etc.
%   bc      - type of bc 0:u=0, 1: du/dx=iku 
%
% output:
%   H       - sparse matrix
%
% Tristan van Leeuwen, 2013
% tristan.vanleeuwen@gmail.com

%% determine dimension and number of gridpoints
nd        = sum(n>1);
N         = prod(n);

% set size of singular dimensions
n(nd+1:3) = 1;
d(nd+1:3) = d(1);

%% FD coefficients first dimension
c1  = ones(N,1)*[1 -2 1]/d(1)^2;

% set bc's
Ia1       = 1:n(1):N;
Ib1       = n(1):n(1):N;        
c1(Ia1,3) = 0; c1(Ia1,2) = c1(Ia1,2) + bc*(1-1i*k(Ia1(:))*d(1))/d(1)^2;
c1(Ib1,1) = 0; c1(Ib1,2) = c1(Ib1,2) + bc*(1-1i*k(Ib1(:))*d(1))/d(1)^2;

%% FD coefficients second dimension
if nd > 1
    c2  = ones(N,1)*[1 -2 1]/d(2)^2;
    
    % set bc's
    Ia2       = ones(n(3),1)*[1:n(1)] + [0:n(1)*n(2):N-n(1)*n(2)]'*ones(1,n(1));
    Ib2       = ones(n(3),1)*[1:n(1)] + [(n(2)-1)*n(1):n(1)*n(2):N-n(1)]'*ones(1,n(1));
    c2(Ia2,3) = 0; c2(Ia2,2) = c2(Ia2,2) + bc*(1-1i*k(Ia2(:))*d(2))/d(2)^2;
    c2(Ib2,1) = 0; c2(Ib2,2) = c2(Ib2,2) + bc*(1-1i*k(Ib2(:))*d(2))/d(2)^2;
end

%% FD coefficients third dimension
if nd > 2
    c3  = ones(N,1)*[1 -2 1]/d(3)^2;
    
    % set bc's
    Ia3       = 1:n(1)*n(2);
    Ib3       = [1:n(1)*n(2)] + n(1)*n(2)*(n(3)-1);
    c3(Ia3,3) = 0; c3(Ia3,2) = c3(Ia3,2) + bc*(1-1i*k(Ia3(:))*d(3))/d(3)^2;
    c3(Ib3,1) = 0; c3(Ib3,2) = c3(Ib3,2) + bc*(1-1i*k(Ib3(:))*d(3))/d(3)^2;
end


%% construct matrix
H = spdiags(k.^2,0,N,N);
H = H + spdiags(c1,[-1 0 1],N,N);
if nd > 1
    H = H + spdiags(c2,[-n(1) 0 n(1)],N,N);
end
if nd > 2
    H = H + spdiags(c3,[-n(1)*n(2) 0 n(1)*n(2)],N,N);
end