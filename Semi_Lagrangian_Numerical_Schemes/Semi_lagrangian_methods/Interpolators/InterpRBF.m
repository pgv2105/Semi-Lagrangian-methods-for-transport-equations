
function W = InterpRBF(X,X0,mPHS,pd)

warning('off')

xst = X(:,1)-X0(1);
yst = X(:,2)-X0(2);
n = length(xst); 

% distance b/w nodes in stencil:
[xx, xc] = meshgrid(X(:,1));
[yy, yc] = meshgrid(X(:,2));
rx = (xx-xc);
ry = (yy-yc);
r2 = rx.^2 + ry.^2;

% eval. nodes:
r2e = xst.^2 + yst.^2;

% RBF-FD collocation matrix & rhs
A = r2.^(mPHS/2);
b = r2e.^(mPHS/2);

% Poly rhs
P = Polybasis(pd,xst,yst);
aux = zeros(size(P,2),1);
aux(1) = 1;

% Augmented RBF-FD collocation matrix & rhs
A = [A P; P' zeros(size(P,2))];
f = [b; aux];

%RBF-FD weights
w = A\f;
W = w(1:n,1);

end




%-----------------------------------------------------------------------
% Polynomial basis for computation of RBF-FD weights

% function [P,fx,fy,f3x2yx,f2xy3y,flap2,flap,f2x,fxy]=Polybasis(n,xst,yst)
function [P]=Polybasis(n,xst,yst)



if n == -1
    
     P = [];
    fx = [];
    fy = [];
f3x2yx = [];
f2xy3y = [];
 flap2 = [];
  flap = [];
   f2x = [];
   fxy = [];

elseif n==0
    
     P = ones(size(xst));
    fx = zeros(size(P,2),1);
    fy = zeros(size(P,2),1);
f3x2yx = zeros(size(P,2),1);
f2xy3y = zeros(size(P,2),1);
 flap2 = zeros(size(P,2),1);
  flap = zeros(size(P,2),1);
   f2x = zeros(size(P,2),1);
   fxy = zeros(size(P,2),1);

elseif n==1
    
     P = [ones(size(xst))...
          xst yst];
      
    fx = [0; 1; 0];
    fy = [0; 0; 1];
f3x2yx = zeros(size(P,2),1);
f2xy3y = zeros(size(P,2),1);
 flap2 = zeros(size(P,2),1);
  flap = zeros(size(P,2),1); 
   f2x = zeros(size(P,2),1);
   fxy = zeros(size(P,2),1); 

elseif n==2
    
     P = [ones(size(xst))...
          xst yst... 
          xst.^2 xst.*yst yst.^2];
      
    fx = [0; 1; 0; zeros(size(P,2)-3,1)];
    fy = [0; 0; 1; zeros(size(P,2)-3,1)];
f3x2yx = zeros(size(P,2),1);
f2xy3y = zeros(size(P,2),1);
 flap2 = zeros(size(P,2),1);
  flap = [0; 0; 0; 2; 0; 2; zeros(size(P,2)-6,1)];
   f2x = [0; 0; 0; 2; 0; 0; zeros(size(P,2)-6,1)];
   fxy = [0; 0; 0; 0; 1; 0; zeros(size(P,2)-6,1)];
 
elseif n==3
    
     P = [ones(size(xst))...
          xst yst...
          xst.^2 xst.*yst yst.^2 ...
          xst.^3 xst.^2.*yst xst.*yst.^2 yst.^3];
    
    fx = [0; 1; 0; zeros(size(P,2)-3,1)];
    fy = [0; 0; 1; zeros(size(P,2)-3,1)];
f3x2yx = [zeros(6,1); 6; 0; 2; 0];
f2xy3y = [zeros(6,1); 0; 2; 0; 6];
 flap2 = zeros(size(P,2),1);
  flap = [0; 0; 0; 2; 0; 2; zeros(size(P,2)-6,1)];
   f2x = [0; 0; 0; 2; 0; 0; zeros(size(P,2)-6,1)];
   fxy = [0; 0; 0; 0; 1; 0; zeros(size(P,2)-6,1)];
   
elseif n==4
    
     P = [ones(size(xst))...
          xst yst...
          xst.^2 xst.*yst yst.^2 ...
          xst.^3 xst.^2.*yst xst.*yst.^2 yst.^3 ...
          xst.^4 xst.^3.*yst xst.^2.*yst.^2 xst.*yst.^3 yst.^4];
    
%     fx = [0; 1; 0; zeros(size(P,2)-3,1)];
%     fy = [0; 0; 1; zeros(size(P,2)-3,1)];
% f3x2yx = [zeros(6,1); 6; 0; 2; 0; zeros(size(P,2)-10,1)];
% f2xy3y = [zeros(6,1); 0; 2; 0; 6; zeros(size(P,2)-10,1)];
%  flap2 = [zeros(10,1); 24; 0; 8; 0; 24];
%   flap = [0; 0; 0; 2; 0; 2; zeros(size(P,2)-6,1)]; 
%    f2x = [0; 0; 0; 2; 0; 0; zeros(size(P,2)-6,1)];
%    fxy = [0; 0; 0; 0; 1; 0; zeros(size(P,2)-6,1)];
   
elseif n==5
    
     P = [ones(size(xst))...
          xst yst... 
          xst.^2 xst.*yst yst.^2 ...
          xst.^3 xst.^2.*yst xst.*yst.^2 yst.^3 ...
          xst.^4 xst.^3.*yst xst.^2.*yst.^2 xst.*yst.^3 yst.^4 ...
          xst.^5 xst.^4.*yst xst.^3.*yst.^2 xst.^2.*yst.^3 xst.*yst.^4 yst.^5];
    
%     fx = [0; 1; 0; zeros(size(P,2)-3,1)];
%     fy = [0; 0; 1; zeros(size(P,2)-3,1)];
% f3x2yx = [zeros(6,1); 6; 0; 2; 0; zeros(size(P,2)-10,1)];
% f2xy3y = [zeros(6,1); 0; 2; 0; 6; zeros(size(P,2)-10,1)];
%  flap2 = [zeros(10,1); 24; 0; 8; 0; 24; zeros(size(P,2)-15,1)];
%   flap = [0; 0; 0; 2; 0; 2; zeros(size(P,2)-6,1)]; 
%    f2x = [0; 0; 0; 2; 0; 0; zeros(size(P,2)-6,1)];
%    fxy = [0; 0; 0; 0; 1; 0; zeros(size(P,2)-6,1)];
   
elseif n==6
    
     P = [ones(size(xst))...
          xst yst...
          xst.^2 xst.*yst yst.^2 ...
          xst.^3 xst.^2.*yst xst.*yst.^2 yst.^3 ...
          xst.^4 xst.^3.*yst xst.^2.*yst.^2 xst.*yst.^3 yst.^4 ...
          xst.^5 xst.^4.*yst xst.^3.*yst.^2 xst.^2.*yst.^3 xst.*yst.^4 yst.^5 ...
          xst.^6 xst.^5.*yst xst.^4.*yst.^2 xst.^3.*yst.^3 xst.^2.*yst.^4 xst.*yst.^5 yst.^6];
    
    fx = [0; 1; 0; zeros(size(P,2)-3,1)];
    fy = [0; 0; 1; zeros(size(P,2)-3,1)];
f3x2yx = [zeros(6,1); 6; 0; 2; 0; zeros(size(P,2)-10,1)];
f2xy3y = [zeros(6,1); 0; 2; 0; 6; zeros(size(P,2)-10,1)];
 flap2 = [zeros(10,1); 24; 0; 8; 0; 24; zeros(size(P,2)-15,1)];
  flap = [0; 0; 0; 2; 0; 2; zeros(size(P,2)-6,1)]; 
   f2x = [0; 0; 0; 2; 0; 0; zeros(size(P,2)-6,1)];
   fxy = [0; 0; 0; 0; 1; 0; zeros(size(P,2)-6,1)];
   
elseif n==7
    
     P = [ones(size(xst))...
          xst yst...
          xst.^2 xst.*yst yst.^2 ...
          xst.^3 xst.^2.*yst xst.*yst.^2 yst.^3 ...
          xst.^4 xst.^3.*yst xst.^2.*yst.^2 xst.*yst.^3 yst.^4 ...
          xst.^5 xst.^4.*yst xst.^3.*yst.^2 xst.^2.*yst.^3 xst.*yst.^4 yst.^5 ...
          xst.^6 xst.^5.*yst xst.^4.*yst.^2 xst.^3.*yst.^3 xst.^2.*yst.^4 xst.*yst.^5 yst.^6 ...
          xst.^7 xst.^6.*yst xst.^5.*yst.^2 xst.^4.*yst.^3 xst.^3.*yst.^4 xst.^2.*yst.^5 xst.*yst.^6 yst.^7];
    
    fx = [0; 1; 0; zeros(size(P,2)-3,1)];
    fy = [0; 0; 1; zeros(size(P,2)-3,1)];
f3x2yx = [zeros(6,1); 6; 0; 2; 0; zeros(size(P,2)-10,1)];
f2xy3y = [zeros(6,1); 0; 2; 0; 6; zeros(size(P,2)-10,1)];
 flap2 = [zeros(10,1); 24; 0; 8; 0; 24; zeros(size(P,2)-15,1)];
  flap = [0; 0; 0; 2; 0; 2; zeros(size(P,2)-6,1)]; 
   f2x = [0; 0; 0; 2; 0; 0; zeros(size(P,2)-6,1)];
   fxy = [0; 0; 0; 0; 1; 0; zeros(size(P,2)-6,1)];
   
elseif n==8
    
     P = [ones(size(xst)) ...
          xst yst ...
          xst.^2 xst.*yst yst.^2 ...
          xst.^3 xst.^2.*yst xst.*yst.^2 yst.^3 ...
          xst.^4 xst.^3.*yst xst.^2.*yst.^2 xst.*yst.^3 yst.^4 ...
          xst.^5 xst.^4.*yst xst.^3.*yst.^2 xst.^2.*yst.^3 xst.*yst.^4 yst.^5 ...
          xst.^6 xst.^5.*yst xst.^4.*yst.^2 xst.^3.*yst.^3 xst.^2.*yst.^4 xst.*yst.^5 yst.^6 ...
          xst.^7 xst.^6.*yst xst.^5.*yst.^2 xst.^4.*yst.^3 xst.^3.*yst.^4 xst.^2.*yst.^5 xst.*yst.^6 yst.^7 ...
          xst.^8 xst.^7.*yst xst.^6.*yst.^2 xst.^5.*yst.^3 xst.^4.*yst.^4 xst.^3.*yst.^5 xst.^2.*yst.^6 xst.*yst.^7 yst.^8];
    
    fx = [0; 1; 0; zeros(size(P,2)-3,1)];
    fy = [0; 0; 1; zeros(size(P,2)-3,1)];
f3x2yx = [zeros(6,1); 6; 0; 2; 0; zeros(size(P,2)-10,1)];
f2xy3y = [zeros(6,1); 0; 2; 0; 6; zeros(size(P,2)-10,1)];
 flap2 = [zeros(10,1); 24; 0; 8; 0; 24; zeros(size(P,2)-15,1)];
  flap = [0; 0; 0; 2; 0; 2; zeros(size(P,2)-6,1)]; 
   f2x = [0; 0; 0; 2; 0; 0; zeros(size(P,2)-6,1)];
   fxy = [0; 0; 0; 0; 1; 0; zeros(size(P,2)-6,1)];   
   
elseif n==9
    
     P = [ones(size(xst)) ...
          xst yst ...
          xst.^2 xst.*yst yst.^2 ...
          xst.^3 xst.^2.*yst xst.*yst.^2 yst.^3 ...
          xst.^4 xst.^3.*yst xst.^2.*yst.^2 xst.*yst.^3 yst.^4 ...
          xst.^5 xst.^4.*yst xst.^3.*yst.^2 xst.^2.*yst.^3 xst.*yst.^4 yst.^5 ...
          xst.^6 xst.^5.*yst xst.^4.*yst.^2 xst.^3.*yst.^3 xst.^2.*yst.^4 xst.*yst.^5 yst.^6 ...
          xst.^7 xst.^6.*yst xst.^5.*yst.^2 xst.^4.*yst.^3 xst.^3.*yst.^4 xst.^2.*yst.^5 xst.*yst.^6 yst.^7 ...
          xst.^8 xst.^7.*yst xst.^6.*yst.^2 xst.^5.*yst.^3 xst.^4.*yst.^4 xst.^3.*yst.^5 xst.^2.*yst.^6 xst.*yst.^7 yst.^8 ...
          xst.^9 xst.^8.*yst xst.^7.*yst.^2 xst.^6.*yst.^3 xst.^5.*yst.^4 xst.^4.*yst.^5 xst.^3.*yst.^6 xst.^2.*yst.^7 xst.*yst.^8 yst.^9];
    
    fx = [0; 1; 0; zeros(size(P,2)-3,1)];
    fy = [0; 0; 1; zeros(size(P,2)-3,1)];
f3x2yx = [zeros(6,1); 6; 0; 2; 0; zeros(size(P,2)-10,1)];
f2xy3y = [zeros(6,1); 0; 2; 0; 6; zeros(size(P,2)-10,1)];
 flap2 = [zeros(10,1); 24; 0; 8; 0; 24; zeros(size(P,2)-15,1)];
  flap = [0; 0; 0; 2; 0; 2; zeros(size(P,2)-6,1)]; 
   f2x = [0; 0; 0; 2; 0; 0; zeros(size(P,2)-6,1)];
   fxy = [0; 0; 0; 0; 1; 0; zeros(size(P,2)-6,1)];
   
elseif n==10
    
     P = [ones(size(xst)) ...
          xst yst ...
          xst.^2 xst.*yst yst.^2 ...
          xst.^3 xst.^2.*yst xst.*yst.^2 yst.^3 ...
          xst.^4 xst.^3.*yst xst.^2.*yst.^2 xst.*yst.^3 yst.^4 ...
          xst.^5 xst.^4.*yst xst.^3.*yst.^2 xst.^2.*yst.^3 xst.*yst.^4 yst.^5 ...
          xst.^6 xst.^5.*yst xst.^4.*yst.^2 xst.^3.*yst.^3 xst.^2.*yst.^4 xst.*yst.^5 yst.^6 ...
          xst.^7 xst.^6.*yst xst.^5.*yst.^2 xst.^4.*yst.^3 xst.^3.*yst.^4 xst.^2.*yst.^5 xst.*yst.^6 yst.^7 ...
          xst.^8 xst.^7.*yst xst.^6.*yst.^2 xst.^5.*yst.^3 xst.^4.*yst.^4 xst.^3.*yst.^5 xst.^2.*yst.^6 xst.*yst.^7 yst.^8 ...
          xst.^9 xst.^8.*yst xst.^7.*yst.^2 xst.^6.*yst.^3 xst.^5.*yst.^4 xst.^4.*yst.^5 xst.^3.*yst.^6 xst.^2.*yst.^7 xst.*yst.^8 yst.^9 ...
          xst.^10 xst.^9.*yst xst.^8.*yst.^2 xst.^7.*yst.^3 xst.^6.*yst.^4 xst.^5.*yst.^5 xst.^4.*yst.^6 xst.^3.*yst.^7 xst.^2.*yst.^8 xst.*yst.^9 yst.^10];
    
    fx = [0; 1; 0; zeros(size(P,2)-3,1)];
    fy = [0; 0; 1; zeros(size(P,2)-3,1)];
f3x2yx = [zeros(6,1); 6; 0; 2; 0; zeros(size(P,2)-10,1)];
f2xy3y = [zeros(6,1); 0; 2; 0; 6; zeros(size(P,2)-10,1)];
 flap2 = [zeros(10,1); 24; 0; 8; 0; 24; zeros(size(P,2)-15,1)];
  flap = [0; 0; 0; 2; 0; 2; zeros(size(P,2)-6,1)]; 
   f2x = [0; 0; 0; 2; 0; 0; zeros(size(P,2)-6,1)];
   fxy = [0; 0; 0; 0; 1; 0; zeros(size(P,2)-6,1)];
   
elseif n==11
    
     P = [ones(size(xst)) ...
          xst yst ...
          xst.^2 xst.*yst yst.^2 ...
          xst.^3 xst.^2.*yst xst.*yst.^2 yst.^3 ...
          xst.^4 xst.^3.*yst xst.^2.*yst.^2 xst.*yst.^3 yst.^4 ...
          xst.^5 xst.^4.*yst xst.^3.*yst.^2 xst.^2.*yst.^3 xst.*yst.^4 yst.^5 ...
          xst.^6 xst.^5.*yst xst.^4.*yst.^2 xst.^3.*yst.^3 xst.^2.*yst.^4 xst.*yst.^5 yst.^6 ...
          xst.^7 xst.^6.*yst xst.^5.*yst.^2 xst.^4.*yst.^3 xst.^3.*yst.^4 xst.^2.*yst.^5 xst.*yst.^6 yst.^7 ...
          xst.^8 xst.^7.*yst xst.^6.*yst.^2 xst.^5.*yst.^3 xst.^4.*yst.^4 xst.^3.*yst.^5 xst.^2.*yst.^6 xst.*yst.^7 yst.^8 ...
          xst.^9 xst.^8.*yst xst.^7.*yst.^2 xst.^6.*yst.^3 xst.^5.*yst.^4 xst.^4.*yst.^5 xst.^3.*yst.^6 xst.^2.*yst.^7 xst.*yst.^8 yst.^9 ...
          xst.^10 xst.^9.*yst xst.^8.*yst.^2 xst.^7.*yst.^3 xst.^6.*yst.^4 xst.^5.*yst.^5 xst.^4.*yst.^6 xst.^3.*yst.^7 xst.^2.*yst.^8 xst.*yst.^9 yst.^10 ...
          xst.^11 xst.^10.*yst xst.^9.*yst.^2 xst.^8.*yst.^3 xst.^7.*yst.^4 xst.^6.*yst.^5 xst.^5.*yst.^6 xst.^4.*yst.^7 xst.^3.*yst.^8 xst.^2.*yst.^9 xst.*yst.^10 yst.^11];
    
    fx = [0; 1; 0; zeros(size(P,2)-3,1)];
    fy = [0; 0; 1; zeros(size(P,2)-3,1)];
f3x2yx = [zeros(6,1); 6; 0; 2; 0; zeros(size(P,2)-10,1)];
f2xy3y = [zeros(6,1); 0; 2; 0; 6; zeros(size(P,2)-10,1)];
 flap2 = [zeros(10,1); 24; 0; 8; 0; 24; zeros(size(P,2)-15,1)];
  flap = [0; 0; 0; 2; 0; 2; zeros(size(P,2)-6,1)]; 
   f2x = [0; 0; 0; 2; 0; 0; zeros(size(P,2)-6,1)];
   fxy = [0; 0; 0; 0; 1; 0; zeros(size(P,2)-6,1)];
   
elseif n==12
    
     P = [ones(size(xst)) ...
          xst yst ...
          xst.^2 xst.*yst yst.^2 ...
          xst.^3 xst.^2.*yst xst.*yst.^2 yst.^3 ...
          xst.^4 xst.^3.*yst xst.^2.*yst.^2 xst.*yst.^3 yst.^4 ...
          xst.^5 xst.^4.*yst xst.^3.*yst.^2 xst.^2.*yst.^3 xst.*yst.^4 yst.^5 ...
          xst.^6 xst.^5.*yst xst.^4.*yst.^2 xst.^3.*yst.^3 xst.^2.*yst.^4 xst.*yst.^5 yst.^6 ...
          xst.^7 xst.^6.*yst xst.^5.*yst.^2 xst.^4.*yst.^3 xst.^3.*yst.^4 xst.^2.*yst.^5 xst.*yst.^6 yst.^7 ...
          xst.^8 xst.^7.*yst xst.^6.*yst.^2 xst.^5.*yst.^3 xst.^4.*yst.^4 xst.^3.*yst.^5 xst.^2.*yst.^6 xst.*yst.^7 yst.^8 ...
          xst.^9 xst.^8.*yst xst.^7.*yst.^2 xst.^6.*yst.^3 xst.^5.*yst.^4 xst.^4.*yst.^5 xst.^3.*yst.^6 xst.^2.*yst.^7 xst.*yst.^8 yst.^9 ...
          xst.^10 xst.^9.*yst xst.^8.*yst.^2 xst.^7.*yst.^3 xst.^6.*yst.^4 xst.^5.*yst.^5 xst.^4.*yst.^6 xst.^3.*yst.^7 xst.^2.*yst.^8 xst.*yst.^9 yst.^10 ...
          xst.^11 xst.^10.*yst xst.^9.*yst.^2 xst.^8.*yst.^3 xst.^7.*yst.^4 xst.^6.*yst.^5 xst.^5.*yst.^6 xst.^4.*yst.^7 xst.^3.*yst.^8 xst.^2.*yst.^9 xst.*yst.^10 yst.^11 ...
          xst.^12 xst.^11.*yst xst.^10.*yst.^2 xst.^9.*yst.^3 xst.^8.*yst.^4 xst.^7.*yst.^5 xst.^6.*yst.^6 xst.^5.*yst.^7 xst.^4.*yst.^8 xst.^3.*yst.^9 xst.^2.*yst.^10 xst.*yst.^11 yst.^12];
    
    fx = [0; 1; 0; zeros(size(P,2)-3,1)];
    fy = [0; 0; 1; zeros(size(P,2)-3,1)];
f3x2yx = [zeros(6,1); 6; 0; 2; 0; zeros(size(P,2)-10,1)];
f2xy3y = [zeros(6,1); 0; 2; 0; 6; zeros(size(P,2)-10,1)];
 flap2 = [zeros(10,1); 24; 0; 8; 0; 24; zeros(size(P,2)-15,1)];
  flap = [0; 0; 0; 2; 0; 2; zeros(size(P,2)-6,1)]; 
   f2x = [0; 0; 0; 2; 0; 0; zeros(size(P,2)-6,1)];
   fxy = [0; 0; 0; 0; 1; 0; zeros(size(P,2)-6,1)];
   
elseif n==13
    
    P = [ones(size(xst)) ...
        xst yst ...
        xst.^2 xst.*yst yst.^2 ...
        xst.^3 xst.^2.*yst xst.*yst.^2 yst.^3 ...
        xst.^4 xst.^3.*yst xst.^2.*yst.^2 xst.*yst.^3 yst.^4 ...
        xst.^5 xst.^4.*yst xst.^3.*yst.^2 xst.^2.*yst.^3 xst.*yst.^4 yst.^5 ...
        xst.^6 xst.^5.*yst xst.^4.*yst.^2 xst.^3.*yst.^3 xst.^2.*yst.^4 xst.*yst.^5 yst.^6 ...
        xst.^7 xst.^6.*yst xst.^5.*yst.^2 xst.^4.*yst.^3 xst.^3.*yst.^4 xst.^2.*yst.^5 xst.*yst.^6 yst.^7 ...
        xst.^8 xst.^7.*yst xst.^6.*yst.^2 xst.^5.*yst.^3 xst.^4.*yst.^4 xst.^3.*yst.^5 xst.^2.*yst.^6 xst.*yst.^7 yst.^8 ...
        xst.^9 xst.^8.*yst xst.^7.*yst.^2 xst.^6.*yst.^3 xst.^5.*yst.^4 xst.^4.*yst.^5 xst.^3.*yst.^6 xst.^2.*yst.^7 xst.*yst.^8 yst.^9 ...
        xst.^10 xst.^9.*yst xst.^8.*yst.^2 xst.^7.*yst.^3 xst.^6.*yst.^4 xst.^5.*yst.^5 xst.^4.*yst.^6 xst.^3.*yst.^7 xst.^2.*yst.^8 xst.*yst.^9 yst.^10 ...
        xst.^11 xst.^10.*yst xst.^9.*yst.^2 xst.^8.*yst.^3 xst.^7.*yst.^4 xst.^6.*yst.^5 xst.^5.*yst.^6 xst.^4.*yst.^7 xst.^3.*yst.^8 xst.^2.*yst.^9 xst.*yst.^10 yst.^11 ...
        xst.^12 xst.^11.*yst xst.^10.*yst.^2 xst.^9.*yst.^3 xst.^8.*yst.^4 xst.^7.*yst.^5 xst.^6.*yst.^6 xst.^5.*yst.^7 xst.^4.*yst.^8 xst.^3.*yst.^9 xst.^2.*yst.^10 xst.*yst.^11 yst.^12 ...
        xst.^13 xst.^12.*yst xst.^11.*yst.^2 xst.^10.*yst.^3 xst.^9.*yst.^4 xst.^8.*yst.^5 xst.^7.*yst.^6 xst.^6.*yst.^7 xst.^5.*yst.^8 xst.^4.*yst.^9 xst.^3.*yst.^10 xst.^2.*yst.^11 xst.*yst.^12 yst.^13];
    
    fx = [0; 1; 0; zeros(size(P,2)-3,1)];
    fy = [0; 0; 1; zeros(size(P,2)-3,1)];
f3x2yx = [zeros(6,1); 6; 0; 2; 0; zeros(size(P,2)-10,1)];
f2xy3y = [zeros(6,1); 0; 2; 0; 6; zeros(size(P,2)-10,1)];
 flap2 = [zeros(10,1); 24; 0; 8; 0; 24; zeros(size(P,2)-15,1)];
  flap = [0; 0; 0; 2; 0; 2; zeros(size(P,2)-6,1)]; 
   f2x = [0; 0; 0; 2; 0; 0; zeros(size(P,2)-6,1)];
   fxy = [0; 0; 0; 0; 1; 0; zeros(size(P,2)-6,1)];
   
else
    
    error('The maximum polynomial degree cannot be higher than 13')
 
end

end