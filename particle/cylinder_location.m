% MATLAB code to reproduce Fig. 9
% Yidan Xue, Sep 2023, last update Jun 2024

warning off
% parameters
r = 0.1;
X = linspace(-1,1,41);
Y = linspace(-1,1,41);
g0c_err = zeros(length(Y),length(X));
g1c_err = zeros(length(Y),length(X));
g2c_err = zeros(length(Y),length(X));
g0c_err2 = zeros(length(Y),length(X));
g1c_err2 = zeros(length(Y),length(X));
g2c_err2 = zeros(length(Y),length(X));
% diameters
Dp = 1;

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
P0 = 0;

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
Z0 = [w1,w2,w3,w4,0+0i].';
Z1 = [w7,w8,w9,w1,0+0i].';
Z2 = [w4,w5,w6,w7,0+0i].';

inpoly = @(z,w) inpolygon(real(z),imag(z),real(w),imag(w));
polyareac = @(z) polyarea(real(z),imag(z));

%%
for ii=1:length(X)
    ii/length(X)
    for jj=1:length(Y)
        ctr = X(ii)+1i*Y(jj);
        np = 100;
        Zp = ctr+r*exp(2i*pi*(1:np)'/np);
        j = inpoly(Zp,Zb);
        if ismember(0, j)
            g0c_err(jj,ii) = NaN;
            g1c_err(jj,ii) = NaN;
            g2c_err(jj,ii) = NaN;
            g0c_err2(jj,ii) = NaN;
            g1c_err2(jj,ii) = NaN;
            g2c_err2(jj,ii) = NaN;
        else
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
            npole = 24; nl = 20; n = 80;
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
            P11 = 0; P21 = 1;
            P12 = 1; P22 = 0;
            % first solution
            P1 = P11; P2 = P21;
            [A1,rhs1,A2,rhs2,PSI,U,V,P] = makerows(Z,Hes,ctr,Pol);  
            
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
            A1(lp,:) =   U(lp,:); rhs1(lp) = 0;  
            A2(lp,:) =   V(lp,:); rhs2(lp) = 0;
            A = [A1; A2]; rhs = [rhs1; rhs2];
            
            % solution and plot
            [A,rhs] = rowweighting(A,rhs,Z,w);
            c = A\rhs;                                  % solve least-squares problem
            [psi,uv,p,omega,f,g] = makefuns(c,Hes,ctr,Pol);   % make function handles
    
            Q11 = psi((w9+w1)/2)-psi((w8+w7)/2);
            Q21 = psi((w6+w7)/2)-psi((w5+w4)/2);
            
            % second solution
            P1 = P12; P2 = P22;
            [A1,rhs1,A2,rhs2,PSI,U,V,P] = makerows(Z,Hes,ctr,Pol);  
            
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
            A1(lp,:) =   U(lp,:); rhs1(lp) = 0;  
            A2(lp,:) =   V(lp,:); rhs2(lp) = 0;
            A = [A1; A2]; rhs = [rhs1; rhs2];
            
            % solution and plot
            [A,rhs] = rowweighting(A,rhs,Z,w);
            c = A\rhs;                                  % solve least-squares problem
            [psi,uv,p,omega,f,g] = makefuns(c,Hes,ctr,Pol);   % make function handles
            
            Q12 = psi((w9+w1)/2)-psi((w8+w7)/2);
            Q22 = psi((w6+w7)/2)-psi((w5+w4)/2);
            [G0c,G1c,G2c] = flow_network(P11,P21,P12,P22,Q11,Q12,Q21,Q22);
            G0ci = 1/(12*l); G1ci = D1^3/(12*l); G2ci = D2^3/(12*l);
            g0c_err(jj,ii) = (G0c-G0ci)/G0ci;
            g1c_err(jj,ii) = (G1c-G1ci)/G1ci;
            g2c_err(jj,ii) = (G2c-G2ci)/G2ci;
            % scale the diameter based on available area
            j0 = inpoly(Zp,Z0);
            j1 = inpoly(Zp,Z1);
            j2 = inpoly(Zp,Z2);
            jorigin = inpoly(Zp,0+0i);
            if jorigin == 1
                Zj0 = flip(Zp(j0).');
                [~, sortj0] = sort(mod(angle(Zj0-ctr),2*pi),'descend');
                Zj0 = Zj0(sortj0);
                A0_new = polyareac([w1,w2,w3,w4,Zj0]);
            elseif ismember(1, j0)
                A0_new = polyareac([w1,w2,w3,w4,0+0i])-polyareac(Zp(j0));
            else
                A0_new = polyareac([w1,w2,w3,w4,0+0i]);
            end

            if jorigin == 1
                Zj2 = flip(Zp(j2).');
                [~, sortj2] = sort(mod(angle(Zj2-ctr)-2*pi/3,2*pi),'descend');
                Zj2 = Zj2(sortj2);
                A2_new = polyareac([w4,w5,w6,w7,Zj2]);
            elseif ismember(1, j2)
                A2_new = polyareac([w4,w5,w6,w7,0+0i])-polyareac(Zp(j2));
            else
                A2_new = polyareac([w4,w5,w6,w7,0+0i]);
            end

            if jorigin == 1
                Zj1 = flip(Zp(j1).');
                [~, sortj1] = sort(mod(angle(Zj1-ctr)+2*pi/3,2*pi),'descend');
                Zj1 = Zj1(sortj1);
                A1_new = polyareac([w7,w8,w9,w1,Zj1]);
            elseif ismember(1, j1)
                A1_new = polyareac([w7,w8,w9,w1,0+0i])-polyareac(Zp(j1));
            else
                A1_new = polyareac([w7,w8,w9,w1,0+0i]);
            end
            D0_new = A0_new/l;
            D1_new = A1_new/l;
            D2_new = A2_new/l;
            G0ci2 = D0_new^3/(12*l); G1ci2 = D1_new^3/(12*l); G2ci2 = D2_new^3/(12*l);
            g0c_err2(jj,ii) = (G0c-G0ci2)/G0ci2;
            g1c_err2(jj,ii) = (G1c-G1ci2)/G1ci2;
            g2c_err2(jj,ii) = (G2c-G2ci2)/G2ci2;
        end
    end
end
%%

save('g0_err_r01.mat','g0c_err')
save('g1_err_r01.mat','g1c_err')
save('g2_err_r01.mat','g2c_err')
save('g0_err2_r01.mat','g0c_err2')
save('g1_err2_r01.mat','g1c_err2')
save('g2_err2_r01.mat','g2c_err2')

%%
FS = 'fontsize'; FW = 'fontweight'; NO = 'normal'; LW = 'linewidth';
tiledlayout(2,3,'TileSpacing','tight');
fs = 24;

nexttile
% subplot(1,3,1)
pcolor(X,Y,g0c_err), hold on, colormap(cool)
shading interp, colorbar
plot(Zb([1:end 1]),'k',LW,.8)
title('$(G_{0}-\tilde{G}_{0})/\tilde{G}_{0}$','interpreter','latex',FS,fs)
xlabel('$X_0$','interpreter','latex', FS,fs)
ylabel('$Y_0$','interpreter','latex', FS,fs)
xticks([-1 -0.5 0 0.5 1])
yticks([-1 -0.5 0 0.5 1])
fontsize(gca,fs,'points')
axis square

nexttile
% subplot(1,3,2)
pcolor(X,Y,g1c_err), hold on, colormap(cool)
shading interp, colorbar
plot(Zb([1:end 1]),'k',LW,.8)
title('$(G_{1}-\tilde{G}_{1})/\tilde{G}_{1}$','interpreter','latex',FS,fs)
xlabel('$X_0$','interpreter','latex', FS,fs)
ylabel('$Y_0$','interpreter','latex', FS,fs)
xticks([-1 -0.5 0 0.5 1])
yticks([-1 -0.5 0 0.5 1])
fontsize(gca,fs,'points')
axis square

nexttile
% subplot(1,3,3)
pcolor(X,Y,g2c_err), hold on, colormap(cool)
shading interp, colorbar
plot(Zb([1:end 1]),'k',LW,.8)
title('$(G_{2}-\tilde{G}_{2})/\tilde{G}_{2}$','interpreter','latex',FS,fs)
xlabel('$X_0$','interpreter','latex', FS,fs)
ylabel('$Y_0$','interpreter','latex', FS,fs)
xticks([-1 -0.5 0 0.5 1])
yticks([-1 -0.5 0 0.5 1])
fontsize(gca,fs,'points')
axis square

nexttile
pcolor(X,Y,g0c_err2), hold on, colormap(cool)
shading interp, colorbar
plot(Zb([1:end 1]),'k',LW,.8)
title('$(G_{0}-\hat{G}_{0})/\hat{G}_{0}$','interpreter','latex',FS,fs)
xlabel('$X_0$','interpreter','latex', FS,fs)
ylabel('$Y_0$','interpreter','latex', FS,fs)
xticks([-1 -0.5 0 0.5 1])
yticks([-1 -0.5 0 0.5 1])
fontsize(gca,fs,'points')
axis square

nexttile
pcolor(X,Y,g1c_err2), hold on, colormap(cool)
shading interp, colorbar
plot(Zb([1:end 1]),'k',LW,.8)
title('$(G_{1}-\hat{G}_{1})/\hat{G}_{1}$','interpreter','latex',FS,fs)
xlabel('$X_0$','interpreter','latex', FS,fs)
ylabel('$Y_0$','interpreter','latex', FS,fs)
xticks([-1 -0.5 0 0.5 1])
yticks([-1 -0.5 0 0.5 1])
fontsize(gca,fs,'points')
axis square

nexttile
pcolor(X,Y,g2c_err2), hold on, colormap(cool)
shading interp, colorbar
plot(Zb([1:end 1]),'k',LW,.8)
title('$(G_{2}-\hat{G}_{2})/\hat{G}_{2}$','interpreter','latex',FS,fs)
xlabel('$X_0$','interpreter','latex', FS,fs)
ylabel('$Y_0$','interpreter','latex', FS,fs)
xticks([-1 -0.5 0 0.5 1])
yticks([-1 -0.5 0 0.5 1])
fontsize(gca,fs,'points')
axis square

set(gcf,'units','inches','position',[0,0,18,12])
exportgraphics(gcf,'particle_location_r01.pdf','Resolution',600)

function [G0c,G1c,G2c] = flow_network(P11,P21,P12,P22,Q11,Q12,Q21,Q22)
    A = [-Q11 Q21; -Q12 Q22];
    b = [P11-P21; P12-P22];
    x = A\b;                    % [1/G1c; 1/G2c]
    G1c = 1/x(1);
    G2c = 1/x(2);
    G0c = (-1/((Q11/G1c+P11)/(Q11+Q21))-1/((Q12/G1c+P12)/(Q12+Q22)))/2;
end
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