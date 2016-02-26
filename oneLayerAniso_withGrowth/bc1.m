%boundary conditions
function bc1=F(y1,y2,y3,y4,R,r,n,mu,k1,k2,gamma,g1,g2)
%%%%%%%%%this is from Moulton's paper, if no anisotropic part, this is
%%%%%%%%%right, since every term is multiplied by r^2*R^2;
%C0 = (n^2-1)*(2*R^2*r^2-R^4);
%C1 = 2*R^2*r^3-2*n^2*R^4*r-R^4*r-n^2*r^5;
%C2 = 2*R^2*r^4+2*R^4*r^2;
%C3 = R^4*r^3;

%%%%%%%%%this is derived from Axxxx, which is consistent with the later
%%%%%%%%%anisotropic part from Bxxxx.
%C0 = (n^2-1)*(-R^2/r^2 + 2);% R^2*r^2
%C1 = -R^2/r + 2*r - n^2*(r^3/R^2+2*R^2/r);% R^2*r^2
%C2 = (2*r^2 +2*R^2);%R^2*r^2
%C3 = R^2*r; %R^2*r^2

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

C0 = (n^2-1)*(r*A1212_p+A1212);
C1 = (r*A1212_p-(n^2-1)*A1212 + n^2*(2*A2112+2*A2211-A2222-A1111))*r;
C2 = (r*A1212_p+4*A1212)*r^2;
C3 = A1212*r^3;

%anisotropic part
I4    =(sin(gamma))^2*alpha^-2 + (cos(gamma))^2*alpha^2;
I4_p  =(sin(gamma))^2*EE       + (cos(gamma))^2*FF;
I4_pp =(sin(gamma))^2*EE_p     + (cos(gamma))^2*FF_p; 
W4    =  k1*exp(k2*(-1+I4)^2)*(-1+I4);
W4_p  =  k1*exp(k2*(-1+I4)^2)*I4_p  + W4*2*k2*(-1+I4)*I4_p;
W4_pp =  k1*exp(k2*(-1+I4)^2)*I4_pp + k1*exp(k2*(-1+I4)^2)*I4_p*2*k2*(-1+I4)*I4_p... 
        +W4*2*k2*(-1+I4)*I4_pp      + W4*2*k2*I4_p*I4_p + W4_p*2*k2*(-1+I4)*I4_p;
    
W44    =  k1*exp(k2*(-1+I4)^2)*(1+2*k2*(-1+I4)^2);
W44_p  =  k1*exp(k2*(-1+I4)^2)*4*k2*(-1+I4)*I4_p + W44*2*k2*(-1+I4)*I4_p;
W44_pp =  k1*exp(k2*(-1+I4)^2)*4*k2*(-1+I4)*I4_pp + k1*exp(k2*(-1+I4)^2)*4*k2*I4_p*I4_p + k1*exp(k2*(-1+I4)^2)*4*k2*(-1+I4)*I4_p*2*k2*(-1+I4)*I4_p+...
          W44*2*k2*(-1+I4)*I4_pp + W44*2*k2*I4_p*I4_p + W44_p*2*k2*(-1+I4)*I4_p;

B1111=4*W4*alpha^(-2)*sin(gamma)^2 + 8*W44*alpha^(-4)*sin(gamma)^4;
B1122=8*W44*sin(gamma)^2*cos(gamma)^2;
B1212=4*W4*alpha^(-2)*sin(gamma)^2 + 8*W44*sin(gamma)^2*cos(gamma)^2;
B1221=8*W44*sin(gamma)^2*cos(gamma)^2;
B2112=8*W44*sin(gamma)^2*cos(gamma)^2;
B2121=4*W4*alpha^2*cos(gamma)^2  + 8*W44*cos(gamma)^2*sin(gamma)^2;
B2211=8*W44*sin(gamma)^2*cos(gamma)^2;
B2222=4*W4*alpha^2*cos(gamma)^2  + 8*W44*alpha^4*cos(gamma)^4;

EEEE_p = (8*r^2-20*R^2)*(r^2-R^2)/r^6;
FFFF_p = 12*r^2/R^8*(R^2-2*r^2)*(R^2-r^2);

B1212_p = 4*W4_p*alpha^(-2)*sin(gamma)^2 + 4*W4*EE*sin(gamma)^2 + 8*W44_p*sin(gamma)^2*cos(gamma)^2;
B1212_pp = 4*W4_pp*alpha^(-2)*sin(gamma)^2 + 4*W4_p*EE*sin(gamma)^2 + 4*W4*EE_p*sin(gamma)^2+...
          8*W44_pp*sin(gamma)^2*cos(gamma)^2;
B2112_p = 8*W44_p*sin(gamma)^2*cos(gamma)^2;
B2211_p = 8*W44_p*sin(gamma)^2*cos(gamma)^2;
B1111_p = 4*W4_p*alpha^(-2)*sin(gamma)^2 + 4*W4*EE*sin(gamma)^2 + 8*W44_p*alpha^(-4)*sin(gamma)^4 + 8*W44*EEEE*sin(gamma)^4;
B2222_p = 4*W4_p*alpha^2*cos(gamma)^2    + 4*W4*FF*sin(gamma)^2 + 8*W44_p*alpha^4*cos(gamma)^4    + 8*W44*FFFF*cos(gamma)^4;

C0 = mu*C0 + (n^2-1)*(r*B1212_p+B1212);
C1 = mu*C1 + (r*B1212_p-(n^2-1)*B1212 + n^2*(2*B2112+2*B2211-B2222-B1111))*r;
C2 = mu*C2 + (r*B1212_p+4*B1212)*r^2;
C3 = mu*C3 + B1212*r^3;

bc1 = C3 * y4 + C2 * y3 + C1 * y2 + C0 * y1;
