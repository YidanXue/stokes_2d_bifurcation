% MATLAB code to reproduce Fig. 6
% Yidan Xue, Jul 2023, last update Jun 2024

tic
% parameters
% diameters
Dp = 1;

% angles
alpha = pi/4;
beta = pi/3;
sina = sin(alpha);
cosa = cos(alpha);
tana = tan(alpha);
sinb = sin(beta);
cosb = cos(beta);
tanb = tan(beta);

% inlet length
l = 2;

% variables
D1 = 0.9;
D2 = 0.8;
P0 = 1;
P1 = 0;
P2 = 0;

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
m = 300; s = linspace(-1,1,m);   % clustered pts in (-1,1)
ms = linspace(0,1,m*2); ms1 = ms(1:end/2); ms2 = ms(end/2+1:end);
Z = [(1-ms2).^3*w9+(3*(1-ms2).^2.*ms2+3*(1-ms2).*ms2.^2)*w1+ms2.^3*w2...
    (w2+w3)/2+(w3-w2)/2*s (1-ms).^3*w3+(3*(1-ms).^2.*ms+3*(1-ms).*ms.^2)*w4+ms.^3*w5...
    (w5+w6)/2+(w6-w5)/2*s (1-ms).^3*w6+(3*(1-ms).^2.*ms+3*(1-ms).*ms.^2)*w7+ms.^3*w8...
    (w8+w9)/2+(w9-w8)/2*s (1-ms1).^3*w9+(3*(1-ms1).^2.*ms1+3*(1-ms1).*ms1.^2)*w1+ms1.^3*w2].';   % boundary pts

% indices
l1 = 1:m; 
l2 = m+1:2*m; 
l3 = 2*m+1:3*m;   
l4 = 3*m+1:4*m;
l5 = 4*m+1:5*m;
l6 = 5*m+1:6*m;
l7 = 6*m+1:7*m;
l8 = 7*m+1:8*m;
l9 = 8*m+1:9*m;

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
n = 48;
F = @(z) conj(z);
[r,pol] = aaa(F,Z([l9 l1]),'tol',1e-8);
inpoly = @(z,w) inpolygon(real(z),imag(z),real(w),imag(w));
j = inpoly(pol,Z);
pol1 = pol(~j & angle(w9-w1)<angle(pol-w1)).';
[r,pol] = aaa(F,Z([l3 l4]),'tol',1e-8);
inpoly = @(z,w) inpolygon(real(z),imag(z),real(w),imag(w));
j = inpoly(pol,Z);
pol2 = pol(~j & angle(pol-w4)<angle(w5-w4)).';
[r,pol] = aaa(F,Z([l6 l7]),'tol',1e-8);
inpoly = @(z,w) inpolygon(real(z),imag(z),real(w),imag(w));
j = inpoly(pol,Z);
pol3 = pol(~j & angle(w6-w7)<angle(pol-w7) & angle(pol-w7)<angle(w8-w7)).';
Pol={pol1,pol2,pol3};
Hes = VAorthog(Z,n,Pol);   % Arnoldi Hessenberg matrices

% boundary conditions
[A1,rhs1,A2,rhs2,PSI,U,V,P] = makerows(Z,n,Hes,Pol);

A1(l1,:) =   U(l1,:); rhs1(l1) = 0;  
A2(l1,:) =   V(l1,:); rhs2(l1) = 0;
% A1(l2,:) =   U(l2,:); rhs1(l2) = -6*imag(Z(l2)).^2+1.5;  
A1(l2,:) =   P(l2,:); rhs1(l2) = P0;
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
A = [A1; A2]; rhs = [rhs1; rhs2];

% solution and plot
c = A\rhs;
toc
error = A*c-rhs;
semilogy(abs(error),'.')
grid on, shg

[psi,uv,p,omega,f,g] = makefuns(c,Hes,Pol);
plotcontours(w,Z,psi,uv,Pol)

function [Hes,R] = VAorthog(Z,n,varargin)  % Vand.+Arnoldi orthogonalization
M = length(Z); Pol = []; if nargin == 3, Pol = varargin{1}; end
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
end

function [R0,R1] = VAeval(Z,Hes,varargin)  % Vand.+Arnoldi basis construction
M = length(Z); Pol = []; if nargin == 3, Pol = varargin{1}; end

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
end

function [A1,rhs1,A2,rhs2,PSI,U,V,P] = makerows(Z,n,Hes,varargin)
Pol = []; if nargin == 4, Pol = varargin{1}; end
[R0,R1] = VAeval(Z,Hes,Pol);
M = length(Z); N = 4*size(R0,2); zero = 0*R0;
cZ = spdiags(conj(Z),0,M,M);                    % conj(Z)
PSI = [cZ*R0 R0]; PSI = [imag(PSI) real(PSI)];  % eq (2.1): stream function
U = [cZ*R1-R0 R1]; U = [real(U) -imag(U)];      % eq (3.9): horizontal velocity
V = [-cZ*R1-R0 -R1]; V = [imag(V) real(V)];     % eq (3.9): vertical velocity
P = [4*R1 zero]; P = [real(P) -imag(P)];
A1 = zeros(M,N); rhs1 = zeros(M,1);
A2 = zeros(M,N); rhs2 = zeros(M,1);
end

function [psi,uv,p,omega,f,g] = makefuns(c,Hes,varargin)  % make function handles
Pol = []; if nargin == 3, Pol = varargin{1}; end
cc = c(1:end/2) + 1i*c(end/2+1:end);
reshaper = @(str) @(z) reshape(fh(str,z(:),cc,Hes,Pol),size(z));
  psi = reshaper('psi');    uv = reshaper('uv');    p = reshaper('p');
omega = reshaper('omega');   f = reshaper('f');   g = reshaper('g');
end

function fh = fh(i,Z,cc,Hes,Pol)
[R0,R1] = VAeval(Z,Hes,Pol);
N = size(R0,2);
cf = cc(1:N); cg = cc(N+(1:N));
switch i
   case   'f'  , fh = R0*cf;
   case   'g'  , fh = R0*cg;
   case  'psi' , fh = imag(conj(Z).*(R0*cf) + R0*cg);
   case   'uv' , fh = Z.*conj(R1*cf) - R0*cf + conj(R1*cg); 
   case   'p'  , fh = real(4*R1*cf);                        
   case 'omega', fh = imag(-4*R1*cf);                       
end
end

function plotcontours(w,Z,psi,uv,varargin)   % contour plot
     Pol = []; if nargin == 5, Pol = varargin{1}; end
     MS = 'markersize'; LW = 'linewidth'; CO = 'color';
     x1 = min(real(Z)); x2 = max(real(Z)); xm = mean([x1 x2]); dx = diff([x1 x2]);
     y1 = min(imag(Z)); y2 = max(imag(Z)); ym = mean([y1 y2]); dy = diff([y1 y2]);
     dmax = max(dx,dy); nx = ceil(800*dx/dmax); ny = ceil(800*dy/dmax);
     x = linspace(x1,x2,nx); y = linspace(y1,y2,ny);
     [xx,yy] = meshgrid(x,y); zz = xx + 1i*yy;
     inpolygonc = @(z,w) inpolygon(real(z),imag(z),real(w),imag(w));
     dZw = min(abs(Z-w.'),[],2);  % distance to nearest corner
     ii = find(dZw>5e-3); outside = ~inpolygonc(zz,Z);
     uu = abs(uv(zz)); uu(outside) = NaN; umax = max(max(uu));
     pcolor(x,y,uu), hold on, colormap(gca,parula)
     shading interp, c=colorbar(), caxis([0 umax])
     c.FontSize = 14;  
     c.Label.FontSize = 14;  
     c.Label.String = 'Velocity magnitude';
     plot(Z([1:end 1]),'k',LW,.8)
     pp = psi(zz); pp(outside) = NaN; pmin = min(min(pp)); pmax = max(max(pp));
     lev = pmin+(.1:.1:.9)*(pmax-pmin);
     contour(x,y,pp,lev,'k',LW,.6)
     if nargin==5, plot(cell2mat(Pol),'.r',MS,8), end
     hold off, axis equal off, axis([xm+.7*dx*[-1 1] ym+.7*dy*[-1 1]])
     exportgraphics(gcf,'smooth_example.pdf','Resolution',600)
end
