%% This script reproduces Figure 3 from the paper: 
% T. van Leeuwen - Fourier analysis of the CGMN method for solving the
% Helmholtz equation, ArXiv:1210:2644, 2012.

% define model
v0 = 6000;
v1 = 1000;
n  = 201;
d  = 1/(n-1);

v = v0*ones(n);
v(20:180,20:180) = v0 + (v1-v0)*ceil(phantom('Modified Shepp-Logan',161));

% frequency, 10 gridpoints/wavelength
f = min(v(:))/max(d)/10;

z = d*[0:n-1];
x = z;
[zz,xx] = ndgrid(z,x);

% setup Helmholtz matrix, convert to band-storage format.
A       = HelmND(2*pi*f./v(:),d*[1 1],n*[1 1],1);
[R,idx] = mat2R(A); 
q = zeros(n);q(10,10)=-1/prod(d);q = q(:);
q = exp(-1i*2*pi*f*xx(:)/v0);

a = 1./sqrt(sum(abs(R).^2,2));
R = bsxfun(@times,a,R);
q = bsxfun(@times,a,q);

% define iteration matrices
w1 = 1.5*ones(n^2,1);
w2 = 1*ones(n^2,1);
w3 = polyval(dlmread('p.dat'),v(:)/f/d(1));
w3 = min(w3,1.95);

At1 = @(x)(x - DKSWPR(R,idx,x,zeros(n^2,1),w1));
bt1 = DKSWPR(R,idx,zeros(n^2,1),q,w1);

At2 = @(x)(x - DKSWPR(R,idx,x,zeros(n^2,1),w2));
bt2 = DKSWPR(R,idx,zeros(n^2,1),q,w2);

At3 = @(x)(x - DKSWPR(R,idx,x,zeros(n^2,1),w3));
bt3 = DKSWPR(R,idx,zeros(n^2,1),q,w3);

% run CG
[u1,~,~,niter1,res1] = pcg(At1,bt1,1e-6,1e4);
[u2,~,~,niter2,res2] = pcg(At2,bt2,1e-6,1e4);
[u3,~,~,niter3,res3] = pcg(At3,bt3,1e-6,1e4);

% plot
figure;imagesc(x,z,reshape(2*pi*f./v,n,n));colorbar;
xlabel('x [m]','fontsize',20);ylabel('z [m]','fontsize',20);
set(gca,'fontsize',20);

figure;imagesc(x,z,reshape(real(u1),n,n),.3*[-1 1]*max(abs(u1)));colormap(gray);
xlabel('x [m]','fontsize',20);ylabel('z [m]','fontsize',20);
set(gca,'fontsize',20);

figure;imagesc(x,z,reshape(real(exp(-1i*2*pi*f*xx(:)/v0)),n,n),[-1 1]*max(abs(exp(-1i*2*pi*f*xx(:)/v0))));colormap(gray);
xlabel('x [m]','fontsize',20);ylabel('z [m]','fontsize',20);
set(gca,'fontsize',20);

figure;semilogy(0:niter1,res1(1:niter1+1)/res1(1),0:niter2,res2(1:niter2+1)/res2(1),0:niter3,res3(1:niter3+1)/res3(1),'linewidth',2);
xlabel('iteration','fontsize',20);ylabel('rel. residual','fontsize',20);xlim([0 max(niter1,niter2)]);ylim([1e-6 1])
set(gca,'fontsize',20);

print(1,'-depsc',['../doc/Fig/exp1_k']);
print(2,'-depsc',['../doc/Fig/exp1_u']);
print(3,'-depsc',['../doc/Fig/exp1_q']);
print(4,'-depsc',['../doc/Fig/exp1_r']);
