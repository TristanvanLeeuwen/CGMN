%% This script reproduces Figure 1 from the paper: 
% T. van Leeuwen - Fourier analysis of the CGMN method for solving the
% Helmholtz equation, ArXiv:1210:2644, 2012.

% define amplitude function (cf. eq. (22-26))
% a1
a1 = @(theta,omega,k,h)(-(k.^2-2*h.^(-2)) - 2*h.^(-2).*cos(theta));
% a2
a2 = @(theta,omega,k,h)((k.^2-2*h.^(-2))^2+2*h.^(-4)+omega.*h.^(-2).*(2*(k.^2-2*h.^(-2)).*exp(-1i*theta)+h.^(-2).*exp(-2i*theta)));
% a3
a3 = @(theta,omega,k,h)((k.^2-2*h.^(-2))^2+2*h.^(-4));
% a4
a4 = @(theta,omega,k,h)((k.^2-2*h.^(-2))^2+2*h.^(-4)+omega.*h.^(-2).*(2*(k.^2-2*h.^(-2)).*exp(1i*theta)+h.^(-2).*exp(2i*theta)));
% a
a = @(theta,omega,k,h)(1-omega.*(2-omega).*a1(theta,omega,k,h).^2.*a3(theta,omega,k,h)./(a2(theta,omega,k,h).*a4(theta,omega,k,h)));


% damping parameters
omega = [.01:.01:1.99];

% angles
theta = linspace(-pi,pi,501);

[oo,tt]=ndgrid(omega,theta);

% # of gridpoints per wavelength
ng   = 10;

% evaluation amplitude function for various ng
ng   = 10:.1:30;
wopt = 0*ng;
for l = 1:length(ng)

    Gl = a(tt,oo,1,2*pi/ng(l));
    B  = max(abs(1-Gl),[],2)./min(abs(1-Gl),[],2);
    [~,i] = min(B);
    wopt(l) = omega(i);
end

% evaluation amplitude function for ng = 10;
G10  = a(tt,oo,1,2*pi/10);

% fit polynomial through wopt
p = polyfit(ng,wopt,3);
dlmwrite('p.dat',p);

% plot
figure;imagesc(theta,omega,G10,[0 1]);
xlabel('\theta','fontsize',20);ylabel('\omega','fontsize',20);
set(gca,'fontsize',20);colorbar;

figure;plot(omega,max(abs(1-G10),[],2)./min(abs(1-G10),[],2),'linewidth',2);
xlabel('\omega','fontsize',20);ylabel('max(1-a)/min(1-a)','fontsize',20);
set(gca,'fontsize',20);

figure;plot(ng,wopt,'linewidth',2);
xlabel('n_g','fontsize',20);ylabel('\omega_{opt}','fontsize',20);
set(gca,'fontsize',20);

print(1,'-depsc',['../doc/Fig/G10']);
print(2,'-depsc',['../doc/Fig/kappa']);
print(3,'-depsc',['../doc/Fig/wopt']);




