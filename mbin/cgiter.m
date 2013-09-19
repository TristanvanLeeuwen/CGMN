%% This script reproduces Figure 2 from the paper: 
% T. van Leeuwen - Fourier analysis of the CGMN method for solving the
% Helmholtz equation, ArXiv:1210:2644, 2012.

% damping parameters
omega = [.05:.05:1.95];
% wavenumbers
k     = 50:50:500;

% varying # of gridpoints/wavelength
for l = 1:length(omega)
    for j = 1:length(k)
        % # of gridpoints, gridspacing, etc.
        n = floor(k(j)/(2*pi))*10 + 1;
        h = 1/(n+1);
        ng1(j) = 2*pi/(k(j)*h);
        w = omega(l);
        
        % setup matrix and r.h.s.
        %A  = helm1d(k(j),n,1);
        A  = HelmND(k(j)*ones(n,1),h,n,0);
        b  = zeros(n,1); b(floor(n/2)+1)=1;
        
        % preconditioner
        D  = diag(diag(A*A'));L=tril(A*A',-1);I=speye(n);
        
        At  = @(x)(w*(2-w)*A'*((D+w*L')\(D*((D+w*L)\(A*x)))));
        bt  = w*(2-w)*A'*((D+w*L')\(D*((D+w*L)\b)));
       
        % call CG
        [X,FLAG,RELRES,ITER,RESVEC] = pcg(At,bt,1e-6,1e4);
        niter1(l,j) = ITER;        
    end
end

% fixed # of gridpoints/wavelength
for l = 1:length(omega)
    for j = 1:length(k)
        % # of gridpoints, gridspacing, etc.
        n = floor(max(k)/(2*pi))*10 + 1;
        h = 1/(n+1);
        ng2(j) = 2*pi/(k(j)*h);
        w = omega(l);
        
        % setup matrix and r.h.s.
        A  = helm1d(k(j),n,1);
        b  = zeros(n,1); b(floor(n/2)+1)=1;

        % preconditioner
        D  = diag(diag(A*A'));L=tril(A*A',-1);I=speye(n);
        
        At  = @(x)(w*(2-w)*A'*((D+w*L')\(D*((D+w*L)\(A*x)))));
        bt  = w*(2-w)*A'*((D+w*L')\(D*((D+w*L)\b)));
       
        % call CG
        [X,FLAG,RELRES,ITER] = pcg(At,bt,1e-6,1e4);
        niter2(l,j) = ITER;        
    end
end

% read results from prev experiment.
p = dlmread('p.dat');

% plot
figure;contour(k,omega,niter1,100);hold on;
plot(k,polyval(p,ng1),'k--','linewidth',2);
xlabel('wavenumber [1/m]','fontsize',20);ylabel('\omega','fontsize',20);
set(gca,'fontsize',20);colorbar;


figure;contour(k,omega,niter2,100);hold on;
plot(k,polyval(p,ng2),'k--','linewidth',2);
xlabel('wavenumber [1/m]','fontsize',20);ylabel('\omega','fontsize',20);
set(gca,'fontsize',20);colorbar;

print(1,'-depsc',['../doc/Fig/niter1']);
print(2,'-depsc',['../doc/Fig/niter2']);

