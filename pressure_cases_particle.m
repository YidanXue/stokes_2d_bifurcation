% MATLAB code to reproduce Fig. 12
% Yidan Xue, Jul 2023, last update Apr 2024

Title = {'(a) $P_1-P_2=0$' '(b) $P_1-P_2=10$' '(c) $P_1-P_2=20$'...
    '(d) $P_1-P_2=30$' '(e) $P_1-P_2=35$' '(f) $P_1-P_2=40$'};
p1 = [0 10 20 30 35 40];

MS = 'markersize'; LW = 'linewidth'; CO = 'color'; CB = 'colorbar';
FS = 'fontsize'; FW = 'fontweight'; NO = 'normal'; LW = 'linewidth';
fs = 12;

tiledlayout(3,2,'Padding','tight','TileSpacing','tight');

% parameters
% diameters
Dp = 1;
r = 0.2;

% angles
alpha = pi/4;
beta = pi/4;
sina = sin(alpha);
cosa = cos(alpha);
tana = tan(alpha);
sinb = sin(beta);
cosb = cos(beta);
tanb = tan(beta);

% inlet length
l = 2;

% variables
D1 = 1;
D2 = 1;

% setup - geometry
wc1 = -l;
wc2 = l*cosb-1i*l*sinb;
wc3 = l*cosa+1i*l*sina;
w1 = (Dp/2-D1/(2*cosa))/tana+Dp/2*1i;
w2 = wc1+Dp/2*1i;
w3 = wc1-Dp/2*1i;
w4 = (Dp/2-D2/(2*cosb))/tanb-Dp/2*1i;
w5 = wc2-D2/2*sinb-1i*D2/2*cosb;
w6 = wc2+D2/2*sinb+1i*D2/2*cosb;
w7x = (D1/cosa+D2/cosb)/(2*(tana+tanb));
w7y = w7x*tana-D1/(2*cosa);
w7 = w7x+1i*w7y;
w8 = wc3+D1/2*sina-1i*D1/2*cosa;
w9 = wc3-D1/2*sina+1i*D1/2*cosa;
if alpha == 0
    w1 = Dp/2*1i;
end
if beta == 0
    w4 = -Dp/2*1i;
end
if alpha == pi/2
    w7x = D1/2;
    w7y = -w7x*tanb+D2/(2*cosb);
    w7 = w7x+1i*w7y;
end

w = [w1; w2; w3; w4; w5; w6; w7; w8; w9];   % corners
m = 300; s = tanh(linspace(-14,14,m));   % clustered pts in (-1,1)
Zb = [(w1+w2)/2+(w2-w1)/2*s (w2+w3)/2+(w3-w2)/2*s (w3+w4)/2+(w4-w3)/2*s...
    (w4+w5)/2+(w5-w4)/2*s (w5+w6)/2+(w6-w5)/2*s (w6+w7)/2+(w7-w6)/2*s...
    (w7+w8)/2+(w8-w7)/2*s (w8+w9)/2+(w9-w8)/2*s (w9+w1)/2+(w1-w9)/2*s].';   % boundary pts
np = 100; ctr = 0;
Zp = ctr+r*exp(2i*pi*(1:np)'/np);
Z = [Zp; Zb];

% indices
lp = 1:np;
l1 = 1+np:m+np; 
l2 = m+1+np:2*m+np; 
l3 = 2*m+1+np:3*m+np;   
l4 = 3*m+1+np:4*m+np;
l5 = 4*m+1+np:5*m+np;
l6 = 5*m+1+np:6*m+np;
l7 = 6*m+1+np:7*m+np;
l8 = 7*m+1+np:8*m+np;
l9 = 8*m+1+np:9*m+np;

% anlges of poles
t1 = (angle(w9-w1)+angle(w2-w1))/2;
t2 = (angle(w1-w2)+angle(w3-w2))/2+pi;
t3 = (angle(w2-w3)+angle(w4-w3))/2+pi;
t4 = (angle(w3-w4)+angle(w5-w4))/2+pi;
t5 = (angle(w4-w5)+angle(w6-w5))/2+pi;
t6 = (angle(w5-w6)+angle(w7-w6))/2;
t7 = (angle(w6-w7)+angle(w8-w7))/2;
t8 = (angle(w7-w8)+angle(w9-w8))/2;
t9 = (angle(w8-w9)+angle(w1-w9))/2+pi;

% create poles
nl = 20; npole = 36; n = 80; 
theta = [pi-angle(w9-w1) angle(w5-w4)+pi angle(w8-w7)-angle(w6-w7)];
sigma = zeros(1,3);
for i =1:3
    if theta(i)>pi/2
        sigma(i) = 4;
    else
        sigma(i) = 2;
    end
end
dk1 = l*cluster(npole,sigma(1)); dk2 = l*cluster(npole,sigma(2)); dk3 = l*cluster(npole,sigma(3));
Pol = {w(1)+exp(1i*t1)*dk1,w(4)+exp(1i*t4)*dk2,w(7)+exp(1i*t7)*dk3};   % the poles
Hes = VAorthog(Z,n,ctr,nl,Pol); 

% boundary conditions
[A1,rhs1,A2,rhs2,PSI,U,V,P] = makerows(Z,Hes,ctr,Pol); 

% annotations
ms = linspace(1/4,3/4,m);
Z_rot_1 = [.5*w1+(1-ms).^3*w2+(3*(1-ms).^2.*ms+3*(1-ms).*ms.^2)*w1+ms.^3*w9].';
Z_rot_2 = [.5*w4+(1-ms).^3*w3+(3*(1-ms).^2.*ms+3*(1-ms).*ms.^2)*w4+ms.^3*w5].';
Z_rot_3 = [.5*w7+(1-ms).^3*w8+(3*(1-ms).^2.*ms+3*(1-ms).*ms.^2)*w7+ms.^3*w6].';

for ii = 1:length(p1)
    P1 = p1(ii);
    P2 = 0;

    A1(l1,:) =   U(l1,:); rhs1(l1) = 0;  
    A2(l1,:) =   V(l1,:); rhs2(l1) = 0;
    A1(l2,:) =   U(l2,:); rhs1(l2) = -6*imag(Z(l2)).^2+1.5;  
    % A1(l2,:) =   P(l2,:); rhs1(l2) = P0;
    A2(l2,:) =   V(l2,:); rhs2(l2) = 0;
    A1(l3,:) =   U(l3,:); rhs1(l3) = 0;  
    A2(l3,:) =   V(l3,:); rhs2(l3) = 0; 
    A1(l4,:) =   U(l4,:); rhs1(l4) = 0;  
    A2(l4,:) =   V(l4,:); rhs2(l4) = 0;
    A1(l5,:) =   P(l5,:); rhs1(l5) = P2;
    if tanb>1e16
        A2(l5,:) =   U(l5,:); rhs2(l5) = 0;
    else
        A2(l5,:) =   cosb*V(l5,:)+sinb*U(l5,:); rhs2(l5) = 0;
    end
    A1(l6,:) =   U(l6,:); rhs1(l6) = 0;  
    A2(l6,:) =   V(l6,:); rhs2(l6) = 0;
    A1(l7,:) =   U(l7,:); rhs1(l7) = 0;  
    A2(l7,:) =   V(l7,:); rhs2(l7) = 0;
    A1(l8,:) =   P(l8,:); rhs1(l8) = P1;
    if tana>1e16
        A2(l8,:) =   U(l8,:); rhs2(l8) = 0; 
    else
        A2(l8,:) =   cosa*V(l8,:)-sina*U(l8,:); rhs2(l8) = 0; 
    end
    A1(l9,:) =   U(l9,:); rhs1(l9) = 0;  
    A2(l9,:) =   V(l9,:); rhs2(l9) = 0;
    A1(lp,:) =   U(lp,:); rhs1(lp) = 0;  
    A2(lp,:) =   V(lp,:); rhs2(lp) = 0;
    A = [A1; A2]; rhs = [rhs1; rhs2];
    
    % solution and plot
    [A,rhs] = rowweighting(A,rhs,Z,w);
    c = A\rhs;                                  % solve least-squares problem
    error = A*c-rhs;
    max(abs(error))
    [psi,uv,p,omega,f,g] = makefuns(c,Hes,ctr,Pol);   % make function handles

     x1 = min(real(Z)); x2 = max(real(Z)); xm = mean([x1 x2]); dx = diff([x1 x2]);
     y1 = min(imag(Z)); y2 = max(imag(Z)); ym = mean([y1 y2]); dy = diff([y1 y2]);
     dmax = max(dx,dy); nx = ceil(400*dx/dmax); ny = ceil(400*dy/dmax);
     x = linspace(x1,x2,nx); y = linspace(y1,y2,ny);
    [xx,yy] = meshgrid(x,y); zz = xx + 1i*yy;
    inpolygonc = @(z,w) inpolygon(real(z),imag(z),real(w),imag(w));
    outside = ~inpolygonc(zz,Zb);
    outsidep = inpolygonc(zz,Zp);
    nexttile
    uu = abs(uv(zz)); uu(outside) = NaN; uu(outsidep) = NaN; umax = max(max(uu));
     pcolor(x,y,uu), hold on, colormap(gca,parula)
     shading interp, c=colorbar, caxis([0 umax])
     c.Label.String = 'Velocity magnitude';
     pp = psi(zz); pp(outside) = NaN; pp(outsidep) = NaN; pmin = min(min(pp)); pmax = max(max(pp));
    lev = pmin+(.1:.1:.9)*(pmax-pmin);
     contour(x,y,pp,lev,'k',LW,.6)
     if P1<=33.22
        lev_sep = [psi(w7) psi(w7)];
    else
        lev_sep = [psi(w1) psi(w1)];
    end
     contour(x,y,pp,lev_sep,'r',LW,.8)
     plot(Zb([1:end 1]),'k','linewidth',.8)
    plot(Zp([1:end 1]),'k','linewidth',.8)

    % arrows
    if ii<=4
        plot(real(Z_rot_1),imag(Z_rot_1),'k',LW,1)
        arrow([real(Z_rot_1(end-1)) imag(Z_rot_1(end-1))],[real(Z_rot_1(end)) imag(Z_rot_1(end))],...
        'length',8,'baseangle',45)
    else
        plot(real(Z_rot_3),imag(Z_rot_3),'k',LW,1)
        arrow([real(Z_rot_3(end-1)) imag(Z_rot_3(end-1))],[real(Z_rot_3(end)) imag(Z_rot_3(end))],...
        'length',8,'baseangle',45)
    end
    plot(real(Z_rot_2),imag(Z_rot_2),'k',LW,1)
    arrow([real(Z_rot_2(end-1)) imag(Z_rot_2(end-1))],[real(Z_rot_2(end)) imag(Z_rot_2(end))],...
    'length',8,'baseangle',45)

     title(Title(ii),'interpreter','latex')
     hold off, axis equal off

end
set(gcf,'units','inches','position',[0,0,6,6])
exportgraphics(gcf,'./pressure_cases_particle.pdf','Resolution',600)

function [Hes,R] = VAorthog(Z,n,ctr,nl,varargin)  % Vand.+Arnoldi orthogonalization
% Input:    Z = column vector of sample points
%           n = degree of polynomial (>=0)
%           ctr = centre of the inner cylinder
%           nl = degree of Laurent series
% Output:   Hes = cell array of Hessenberg matrices (length 1+length(Pol))
%           R = matrix of basis vectors
M = length(Z); Pol = []; if nargin == 5, Pol = varargin{1}; end
% First orthogonalize the polynomial part
Q = ones(M,1); H = zeros(n+1,n);
for k = 1:n
    q = Z.*Q(:,k);
    for j = 1:k, H(j,k) = Q(:,j)'*q/M; q = q - H(j,k)*Q(:,j); end
    H(k+1,k) = norm(q)/sqrt(M); Q(:,k+1) = q/H(k+1,k);
end
Hes{1} = H; R = Q;
% Next orthogonalize the pole parts, if any
while ~isempty(Pol)
   pol = Pol{1}; Pol(1) = [];
   np = length(pol); H = zeros(np,np-1); Q = ones(M,1);
   for k = 1:np
      q = Q(:,k)./(Z-pol(k));
      for j = 1:k, H(j,k) = Q(:,j)'*q/M; q = q - H(j,k)*Q(:,j); end
      H(k+1,k) = norm(q)/sqrt(M); Q(:,k+1) = q/H(k+1,k);
    end
   Hes{length(Hes)+1} = H; R = [R Q(:,2:end)];
end
% Next orthogonalize the Laurent series
for m = 1:length(ctr)
    Q = ones(M,1); H = zeros(nl+1,nl);
    for k = 1:nl
       q = 1./(Z-ctr(m)).*Q(:,k);
       for j = 1:k, H(j,k) = Q(:,j)'*q/M; q = q - H(j,k)*Q(:,j); end
       H(k+1,k) = norm(q)/sqrt(M); Q(:,k+1) = q/H(k+1,k);
    end
    Hes{length(Hes)+1} = H; R = [R Q(:,2:end)];
end
end

function [R0,R1] = VAeval(Z,Hes,ctr,varargin)  % Vand.+Arnoldi basis construction
% Input:    Z = column vector of sample points
%           n = degree of polynomial (>=0)
%           ctr = centre of the inner cylinder
% Output:   R0 = matrix of basis vectors for functions
%           R1 = matrix of basis vectors for derivatives
M = length(Z); Pol = []; if nargin == 4, Pol = varargin{1}; end
% First construct the polynomial part of the basis
H = Hes{1}; Hes(1) = []; n = size(H,2);
Q = ones(M,1); D = zeros(M,1);
for k = 1:n
    hkk = H(k+1,k);
    Q(:,k+1) = ( Z.*Q(:,k) - Q(:,1:k)*H(1:k,k))/hkk;
    D(:,k+1) = ( Z.*D(:,k) - D(:,1:k)*H(1:k,k) + Q(:,k) )/hkk;
end
R0 = Q; R1 = D;
% Next construct the pole parts of the basis, if any
while ~isempty(Pol)
    pol = Pol{1}; Pol(1) = [];
    H = Hes{1}; Hes(1) = []; np = length(pol); Q = ones(M,1); D = zeros(M,1);
    for k = 1:np
       Zpki = 1./(Z-pol(k)); hkk = H(k+1,k);
       Q(:,k+1) = ( Q(:,k).*Zpki - Q(:,1:k)*H(1:k,k)                   )/hkk;
       D(:,k+1) = ( D(:,k).*Zpki - D(:,1:k)*H(1:k,k) - Q(:,k).*Zpki.^2 )/hkk;
    end
    R0 = [R0 Q(:,2:end)]; R1 = [R1 D(:,2:end)];
end
% Next construct the basis for the Laurent series
for m = 1:length(ctr)
    H = Hes{1}; Hes(1) = []; nl = size(H,2); Q = ones(M,1); D = zeros(M,1);
    Zpki = 1./(Z-ctr(m)); Zpkid = -1./(Z-ctr(m)).^2;
    for k = 1:nl
        hkk = H(k+1,k);
        Q(:,k+1) = ( Q(:,k).*Zpki - Q(:,1:k)*H(1:k,k))/hkk;
        D(:,k+1) = ( D(:,k).*Zpki - D(:,1:k)*H(1:k,k) + Q(:,k).*Zpkid )/hkk;
    end
    R0 = [R0 Q(:,2:end)]; R1 = [R1 D(:,2:end)];
end
end
function d = cluster(n,sigma)   % n points exponentially clustered in (0,1]
nc = ceil(n); d = exp(sigma*(sqrt(nc:-1:1)-sqrt(nc)));   % originally it was 4
end
function [A,rhs] = rowweighting(A,rhs,Z,w)
dZw = min(abs(Z-w.'),[],2);
wt = [dZw; dZw];
M2 = 2*length(Z); W = spdiags(wt,0,M2,M2);
A = W*A; rhs = W*rhs;
end
function [A1,rhs1,A2,rhs2,PSI,U,V,P] = makerows(Z,Hes,ctr,varargin)
Pol = []; if nargin == 4, Pol = varargin{1}; end
m = length(ctr);
[R0,R1] = VAeval(Z,Hes,ctr,Pol);
M = length(Z); N = 4*size(R0,2)+4*m; zero = 0*R0;
cZ = spdiags(conj(Z),0,M,M);                              % conj(Z)
oZ = 1./(Z-ctr);                                          % 1/(Z-ctr)
lZ = log(Z-ctr);                                          % log(Z-ctr)
PSI = [cZ*R0 R0];                                         % stream function
U = [cZ*R1-R0 R1];                                        % horizontal vel.
V = [-cZ*R1-R0 -R1];                                      % vertical vel.  
P = [4*R1 zero];                                          % pressure
PSI = [imag(PSI) imag(cZ*lZ-(Z-ctr).*lZ+Z) imag(lZ)...
    real(PSI) real(cZ*lZ+(Z-ctr).*lZ-Z) real(lZ)];        % log terms
U = [real(U) real(cZ*oZ-2*lZ) real(oZ)...
    -imag(U) -imag(cZ*oZ) -imag(oZ)];                     % log terms
V = [imag(V) imag(-cZ*oZ) imag(-oZ)...
    real(V) real(-cZ*oZ-2*lZ) real(-oZ)];                 % log terms 
P = [real(P) real(4*oZ) zeros(M,m) ...
    -imag(P) -imag(4*oZ) zeros(M,m)];                     % log terms
A1 = zeros(M,N); rhs1 = zeros(M,1);
A2 = zeros(M,N); rhs2 = zeros(M,1);
end

function [psi,uv,p,omega,f,g] = makefuns(c,Hes,ctr,varargin)  % make function handles
Pol = []; if nargin == 4, Pol = varargin{1}; end
cc = c(1:end/2) + 1i*c(end/2+1:end);
reshaper = @(str) @(z) reshape(fh(str,z(:),cc,Hes,ctr,Pol),size(z));
  psi = reshaper('psi');    uv = reshaper('uv');    p = reshaper('p');
omega = reshaper('omega');   f = reshaper('f');   g = reshaper('g');
end

function fh = fh(i,Z,cc,Hes,ctr,Pol)
m = length(ctr);
[R0,R1] = VAeval(Z,Hes,ctr,Pol);
N = size(R0,2);
cf = cc(1:N); cg = cc(N+(1:N)); clf = cc(2*N+1:2*N+m); clg = cc(2*N+m+1:2*(N+m));
%   clf/clg = coefficient for the logarithmic term in f/g
switch i
   case   'f'  , fh = R0*cf;
       for k=1:m fh = fh+clf(k)*log(Z-ctr(k)); end
   case   'g'  , fh = R0*cg;
       for k=1:m fh = fh+clg(k)*log(Z-ctr(k))-conj(clf(k))*((Z-ctr(k)).*log(Z-ctr(k))-Z); end
   case  'psi' , fh = imag(conj(Z).*(R0*cf)+ R0*cg);
       for k=1:m fh = fh+imag(conj(Z).*(clf(k)*log(Z-ctr(k)))...
           +clg(k)*log(Z-ctr(k))-conj(clf(k))*((Z-ctr(k)).*log(Z-ctr(k))-Z)); end
   case   'uv' , fh = Z.*conj(R1*cf) - R0*cf+ conj(R1*cg);
       for k=1:m fh = fh+Z.*conj(clf(k)./(Z-ctr(k)))- clf(k)*log(Z-ctr(k))...
           + conj(clg(k)./(Z-ctr(k))-conj(clf(k))*log(Z-ctr(k))); end
   case   'p'  , fh = real(4*R1*cf);
       for k=1:m fh = fh+real(4*clf(k)./(Z-ctr(k))); end
   case 'omega', fh = imag(-4*R1*cf); 
       for k=1:m fh = fh+imag(-4*clf(k)./(Z-ctr(k))); end
end
end
