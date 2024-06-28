% MATLAB code to reproduce Fig. 5
% Yidan Xue, Jul 2023, last update Mar 2024

% parameters
num = 41;
Alpha = linspace(0,pi/2,num);
Beta = linspace(0,pi/2,num);
g0c_err = zeros(length(Beta),length(Alpha));
g1c_err = zeros(length(Beta),length(Alpha));
g2c_err = zeros(length(Beta),length(Alpha));
for ii=1:length(Alpha)
    ii/length(Alpha)
    for jj=1:length(Beta)
        % diameters
        Dp = 1;
        
        % angles
        alpha = Alpha(ii);
        beta = Beta(jj);

        if abs(alpha+beta)<pi/2-1e-5
            g0c_err(jj,ii) = NaN;
            g1c_err(jj,ii) = NaN;
            g2c_err(jj,ii) = NaN;
        else

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
            Z = [(w1+w2)/2+(w2-w1)/2*s (w2+w3)/2+(w3-w2)/2*s (w3+w4)/2+(w4-w3)/2*s...
                (w4+w5)/2+(w5-w4)/2*s (w5+w6)/2+(w6-w5)/2*s (w6+w7)/2+(w7-w6)/2*s...
                (w7+w8)/2+(w8-w7)/2*s (w8+w9)/2+(w9-w8)/2*s (w9+w1)/2+(w1-w9)/2*s].';   % boundary pts
            
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
            n = 24; np = 48;
            theta = [pi-angle(w9-w1) angle(w5-w4)+pi angle(w8-w7)-angle(w6-w7)];
            sigma = zeros(1,3);
            for i =1:3
                if theta(i)>pi/2
                    sigma(i) = 4;
                else
                    sigma(i) = 2;
                end
            end
            dk1 = l*cluster(np,sigma(1)); dk2 = l*cluster(np,sigma(2)); dk3 = l*cluster(np,sigma(3));
            Pol = {w(1)+exp(1i*t1)*dk1,w(4)+exp(1i*t4)*dk2,w(7)+exp(1i*t7)*dk3};   % the poles
            
            Hes = VAorthog(Z,n,Pol);   % Arnoldi Hessenberg matrices
            
            % boundary conditions
            P11 = 1; P21 = 0;
            P12 = 0; P22 = 1;
            % first solution
            P1 = P11; P2 = P21;
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
            [A,rhs] = rowweighting(A,rhs,Z,w);
            c = A\rhs;
            
            [psi,uv,p,omega,f,g] = makefuns(c,Hes,Pol);
            % plotcontours(w,Z,psi,uv,Pol)
            Q11 = psi(w9)-psi(w8);
            Q21 = psi(w6)-psi(w5);
            
            % second solution
            P1 = P12; P2 = P22;
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
            [A,rhs] = rowweighting(A,rhs,Z,w);
            c = A\rhs;
            [psi,uv,p,omega,f,g] = makefuns(c,Hes,Pol);
            Q12 = psi(w9)-psi(w8);
            Q22 = psi(w6)-psi(w5);
            
            [G0c,G1c,G2c] = flow_network(P11,P21,P12,P22,Q11,Q12,Q21,Q22);
            G0ci = 1/(12*l); G1ci = D1^3/(12*l); G2ci = D2^3/(12*l);
            g0c_err(jj,ii) = (G0c-G0ci)/G0ci;
            g1c_err(jj,ii) = (G1c-G1ci)/G1ci;
            g2c_err(jj,ii) = (G2c-G2ci)/G2ci;
        end
    end
end
%%
FS = 'fontsize'; FW = 'fontweight'; NO = 'normal'; LW = 'linewidth';
tiledlayout(1,3,'Padding','tight','TileSpacing','tight');
fs = 24;

nexttile
% subplot(1,3,1)
pcolor(Alpha,Beta,g0c_err), hold on, colormap(cool)
emin = min(min(g0c_err));
emax = max(max(g0c_err));
shading interp, c=colorbar(), caxis([0 0.09])
title('$(G_0-\tilde{G}_0)/\tilde{G}_0$','interpreter','latex',FS,fs)
xlabel('$\alpha$','interpreter','latex', FS,fs)
ylabel('$\beta$','interpreter','latex', FS,fs)
xticks([0 1/8 1/4 3/8 1/2]*pi)
yticks([0 1/8 1/4 3/8 1/2]*pi)
set(gca,'XTickLabel',{'0','1/8\pi','1/4\pi','3/8\pi','1/2\pi'})
set(gca,'YTickLabel',{'0','1/8\pi','1/4\pi','3/8\pi','1/2\pi'})
% xlim([pi/4,pi/2])
% ylim([pi/4,pi/2])
fontsize(gca,fs,'points')
axis square

nexttile
% subplot(1,3,2)
pcolor(Alpha,Beta,g1c_err), hold on, colormap(cool)
emin = min(min(g1c_err));
emax = max(max(g1c_err));
shading interp, c=colorbar(), caxis([0 0.09])
title('$(G_1-\tilde{G}_1)/\tilde{G}_1$','interpreter','latex',FS,fs)
xlabel('$\alpha$','interpreter','latex', FS,fs)
ylabel('$\beta$','interpreter','latex', FS,fs)
xticks([0 1/8 1/4 3/8 1/2]*pi)
yticks([0 1/8 1/4 3/8 1/2]*pi)
set(gca,'XTickLabel',{'0','1/8\pi','1/4\pi','3/8\pi','1/2\pi'})
set(gca,'YTickLabel',{'0','1/8\pi','1/4\pi','3/8\pi','1/2\pi'})
% xlim([pi/4,pi/2])
% ylim([pi/4,pi/2])
fontsize(gca,fs,'points')
axis square

nexttile
% subplot(1,3,3)
pcolor(Alpha,Beta,g2c_err), hold on, colormap(cool)
emin = min(min(g2c_err));
emax = max(max(g2c_err));
shading interp, c=colorbar(), caxis([0 0.09])
title('$(G_2-\tilde{G}_2)/\tilde{G}_2$','interpreter','latex',FS,fs)
xlabel('$\alpha$','interpreter','latex', FS,fs)
ylabel('$\beta$','interpreter','latex', FS,fs)
xticks([0 1/8 1/4 3/8 1/2]*pi)
yticks([0 1/8 1/4 3/8 1/2]*pi)
set(gca,'XTickLabel',{'0','1/8\pi','1/4\pi','3/8\pi','1/2\pi'})
set(gca,'YTickLabel',{'0','1/8\pi','1/4\pi','3/8\pi','1/2\pi'})
% xlim([pi/4,pi/2])
% ylim([pi/4,pi/2])
fontsize(gca,fs,'points')
axis square
set(gcf,'units','inches','position',[0,0,18,6])
exportgraphics(gcf,'angle_effects.pdf','Resolution',600)

function [G0c,G1c,G2c] = flow_network(P11,P21,P12,P22,Q11,Q12,Q21,Q22)
    A = [-Q11 Q21; -Q12 Q22];
    b = [P11-P21; P12-P22];
    x = A\b;                    % [1/G1c; 1/G2c]
    G1c = 1/x(1);
    G2c = 1/x(2);
    G0c = (-1/((Q11/G1c+P11)/(Q11+Q21))-1/((Q12/G1c+P12)/(Q12+Q22)))/2;
end

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

function d = cluster(n,sigma)   % n points exponentially clustered in (0,1]
nc = ceil(n); d = exp(sigma*(sqrt(nc:-1:1)-sqrt(nc)));   % originally it was 4
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

function [A,rhs] = rowweighting(A,rhs,Z,w)
dZw = min(abs(Z-w.'),[],2);
wt = [dZw; dZw];
M2 = 2*length(Z); W = spdiags(wt,0,M2,M2);
A = W*A; rhs = W*rhs;
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