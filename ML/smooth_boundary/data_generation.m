%% Create an array with 1000 rows
warning off
para = []; num = 1000;
g0c = []; g1c = []; g2c = [];
L = 2;
while size(para,1)<num
    if mod(size(para,1),100) == 0
        size(para,1)/num
    end
    % generate random parameters
    D1 = 0.5*rand+0.5;
    D2 = 0.5*rand+0.5;
    alpha = rand*pi/2;
    beta = rand*pi/2;

    % parameters
    P0 = 0;
    Dp = 1;

    sina = sin(alpha);
    cosa = cos(alpha);
    tana = tan(alpha);
    sinb = sin(beta);
    cosb = cos(beta);
    tanb = tan(beta);
    
    % setup - geometry
    wc1 = -L;
    wc2 = L*cosb-1i*L*sinb;
    wc3 = L*cosa+1i*L*sina;
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
    % check the geometry is a bifurcation
    if abs(alpha+beta)>=pi/2 && imag(w5)<=imag(w4) && imag(w9)>=imag(w1) && real(w7)<=real(w6) && real(w7)<=real(w8)
        para = [para; D1 D2 alpha beta];
    
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
        c = A\rhs;
        error = A*c-rhs;
        maxerror = norm(error,inf);
        if maxerror>=1e-6
            disp('Error! Parameters:')
            para(end,:)
        end
        
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
        c = A\rhs;
        error = A*c-rhs;
        maxerror = norm(error,inf);
        if maxerror>=1e-6
            disp('Error! Parameters:')
            para(end,:)
        end

        [psi,uv,p,omega,f,g] = makefuns(c,Hes,Pol);
        Q12 = psi(w9)-psi(w8);
        Q22 = psi(w6)-psi(w5);
        
        [G0c,G1c,G2c] = flow_network(P11,P21,P12,P22,Q11,Q12,Q21,Q22);
        g0c = [g0c;G0c];
        g1c = [g1c;G1c];
        g2c = [g2c;G2c];
    end
end

%%
save('parameters.mat','para')
save('g0c.mat','g0c')
save('g1c.mat','g1c')
save('g2c.mat','g2c')

%% functions

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
