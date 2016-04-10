function y4=F(y1,y2,y3,R,r,mu,n,g1,g2,rhs)
%C0 = (n^2-1)*(2*g2^3*R^2*r^2-g1*g2^4*R^4);
%C1 = 2*g2^3*R^2*r^3-2*g1*g2^4*n^2*R^4*r-g1*g2^4*R^4*r-g1*n^2*r^5;
%C2 = (2*g2^3*R^2*r^4+2*g1*g2^4*R^4*r^2);
%C3 = g1*g2^4*R^4*r^3;


PI=3.141592654;
alpha=r/R/g2;
EE = -2*R^2*g2^2/r^3 + 2*g2/r/g1; % d/dr (alpha^(-2))
FF = 2*r/R^2/g2^2 - 2*r^3/R^4/g1/g2^3; % d/dr (alpha^2)
EE_p = -6*g2/r^2/g1 + 6*R^2*g2^2/r^4;
FF_p = 2/R^2/g2^2 - 4*r^2/R^4/g1/g2^3 - 6*r^2*g1/R^4/g2^3 + 8*r^4/R^6/g2^4;
EEEE = 4*R^2/r^3*g2^3/g1 - 4*R^4/r^5*g2^4; % d/dr (alpha^(-4))
FFFF = 4*r^3/R^4/g2^4 - 4*r^5/R^6/g1/g2^5; % d/dr (alpha^4)

%isotropic part
A1111 = alpha^(-2); A1111_p = EE;
A1122 = 0;          
A1212 = alpha^(-2); A1212_p = EE; A1212_pp = EE_p;
A1221 = 0;
A2121 = alpha^2;    
A2112 = 0;          A2112_p = 0;
A2211 = 0;          A2211_p = 0;
A2222 = alpha^2;    A2222_p = FF;

C0 = mu*(n^2-1)*(r*A1212_p+A1212);
C1 = mu*(r*A1212_p-(n^2-1)*A1212 + n^2*(2*A2112+2*A2211-A2222-A1111))*r;
C2 = mu*(r*A1212_p+4*A1212)*r^2;
C3 = mu*A1212*r^3;

y4 = (rhs - (C2 * y3 + C1 * y2 + C0 * y1) )/C3;%*R^2*r^2*g1*g2^2;
