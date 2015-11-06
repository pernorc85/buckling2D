%boundary conditions
function bc1=F(y1,y2,y3,y4,R,r,n,gamma)
C0 = (n^2-1)*(2*R^2*r^2-R^4);
C1 = 2*R^2*r^3-2*n^2*R^4*r-R^4*r-n^2*r^5;
C2 = 2*R^2*r^4+2*R^4*r^2;
C3 = R^4*r^3;

PI=3.141592654;
%gamma=PI/3;%angular orientation of collagen fibers in the cross section
k1=0.1;%parameters for anisotropic material properties
k2=2;
alpha=r/R;
I4    =(sin(gamma))^2*alpha^-2 + (cos(gamma))^2*alpha^2;
I4_p  =(sin(gamma))^2*(-2)*(R^2-r^2)/r^3 + (cos(gamma))^2*(R^2-r^2)*r/R^4;
I4_pp =(sin(gamma))^2*6*(R^2-r^2)/r^4    + (cos(gamma))^2*(R^2-r^2)*(R^2-4*r^2)*2/R^6; 
W4=k1*exp(k2*(-1+I4)^2)*(-1+I4);
W44=k1*exp(k2*(-1+I4)^2)*(1+2*k2*(-1+I4)^2);
W4_p  =  k1*exp(k2*(-1+I4)^2)*I4_p  + W4*2*k2*(-1+I4)*I4_p;
W4_pp =  k1*exp(k2*(-1+I4)^2)*I4_pp + k1*exp(k2*(-1+I4)^2)*I4_p*2*k2*(-1+I4)*I4_p... 
        +W4*2*k2*(-1+I4)*I4_pp      + W4*2*k2*I4_p*I4_p + (k1*exp(k2*(-1+I4)^2)*I4_p + W4*2*k2*(-1+I4)*I4_p) * 2*k2*(-1+I4)*I4_p;
    
W44=k1*exp(k2*(-1+I4)^2)*(1+2*k2*(-1+I4)^2);
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

EE = -2*(R^2 - r^2)/r^3; % d/dr (alpha^(-2))
FF = 2*r*(R^2 - r^2)/R^4; % d/dr (alpha^2)
EEEE = 4*R^2*(r^2 - R^2)/r^5; % d/dr (alpha^(-4))
FFFF = 4*r^3*(R^2 - r^2)/R^6; % d/dr (alpha^4)

EE_p = 6*(R^2 - r^2)/r^4;
FF_p = 2*(R^2 - r^2)*(R^2 - 4*r^2)/R^6;
EEEE_p = (8*r^2-20*R^2)*(r^2-R^2)/r^6;
FFFF_p = 12*r^2/R^8*(R^2-2*r^2)*(R^2-r^2);

B1212_p = 4*W4_p*alpha^(-2)*sin(gamma)^2 + 4*W4*EE*sin(gamma)^2 + 8*W44_p*sin(gamma)^2*cos(gamma)^2;
B1212_pp = 4*W4_pp*alpha^(-2)*sin(gamma)^2 + 4*W4_p*EE*sin(gamma)^2 + 4*W4*EE_p*sin(gamma)^2+...
          8*W44_pp*sin(gamma)^2*cos(gamma)^2;
B2112_p = 8*W44_p*sin(gamma)^2*cos(gamma)^2;
B2211_p = 8*W44_p*sin(gamma)^2*cos(gamma)^2;
B1111_p = 4*W4_p*alpha^(-2)*sin(gamma)^2 + 4*W4*EE*sin(gamma)^2 + 8*W44_p*alpha^(-4)*sin(gamma)^4 + 8*W44*EEEE*sin(gamma)^4;
B2222_p = 4*W4_p*alpha^2*cos(gamma)^2    + 4*W4*FF*sin(gamma)^2 + 8*W44_p*alpha^4*cos(gamma)^4    + 8*W44*FFFF*cos(gamma)^4;

C0 = C0 + (n^2-1)*(r*B1212_p+B1212);
C1 = C1 + (r*B1212_p-(n^2-1)*B1212 + n^2*(2*B2112+2*B2211-B2222-B1111))*r;
C2 = C2 + (r*B1212_p+4*B1212)*r^2;
C3 = C3 + B1212*r^3;

bc1 = C3 * y4 + C2 * y3 + C1 * y2 + C0 * y1;
