%buckling of a isotropic tube without growth
%inner radius after deformation is the only free parameter
function yp=F(r,y,A,a,mu,n,g1,g2)
yp = zeros(5,1);%derivatives

%constants
%n=2;%buckling mode(3 4 5 6 ...)
Rsquare = A^2+(r^2-a^2)/g1/g2;
R       = sqrt(Rsquare);
EE = -2*R^2*g2^2/r^3 + 2*g2/r/g1; % d/dr (alpha^(-2))
FF = 2*r/R^2/g2^2 - 2*r^3/R^4/g1/g2^3; % d/dr (alpha^2)
EE_p = -6*g2/r^2/g1 + 6*R^2*g2^2/r^4;
FF_p = 2/R^2/g2^2 - 4*r^2/R^4/g1/g2^3 - 6*r^2*g1/R^4/g2^3 + 8*r^4/R^6/g2^4;
EEEE = 4*R^2/r^3*g2^3/g1 - 4*R^4/r^5*g2^4; % d/dr (alpha^(-4))
FFFF = 4*r^3/R^4/g2^4 - 4*r^5/R^6/g1/g2^5; % d/dr (alpha^4)

%coefficients
%B0 = (n^2-1)*(3*R^4 - 4*R^2*r^2 + r^4*n^2)/(r^4*R^2);
%B1 = (-3*r^4*n^2*R^2 + 2*r^6*n^2 + 3*R^6 - 4*R^4*r^2 + n^2*R^6 - 2*R^4*n^2*r^2)/(r^3*R^4);
%B2 = -(3*R^4 - 8*R^2*r^2 + n^2*R^4 + r^4*n^2)/(r^2*R^2);
%B3 = 2*(R^2 + 2*r^2)/r;
%B4 = R^2;

alpha=r/R/g2;
I1    = alpha^-2 + alpha^2;
W1    = mu/2;

A1111 = alpha^(-2); A1111_p = EE;
A1122 = 0;          
A1212 = alpha^(-2); A1212_p = EE; A1212_pp = EE_p;
A1221 = 0;
A2121 = alpha^2;    
A2112 = 0;          A2112_p = 0;
A2211 = 0;          A2211_p = 0;
A2222 = alpha^2;    A2222_p = FF;
sigma2=alpha*alpha;       %=alpha*W1*dI1/dalpha;
sigma1=1/alpha*1/alpha;   %=alpha-1*W1*dI1/dalpha-1;
%have to eliminate mu since mu is multiplied later

B0 = (n^2-1)*(r^2*A1212_pp + r*A1212_p + (n^2-1)*A1212 + n^2*(sigma2-sigma1) )/r^2;
%B0 = (n^2-1)*(3*g2^4*g1*(Rsquare)^2 - 4*g2^3*(Rsquare)*r^2 + g1*r^4*n^2)/(r^4*(Rsquare)*g2^2*g1);
B1 = (A1212_p+r*A1212_pp-A1212/r) + n^2/r*(2*A2112+2*A2211-A2222-A1111) + n^2*(2*A2112_p+2*A2211_p-A2222_p-A1111_p);
B2 = (7*r*A1212_p+r^2*A1212_pp+5*A1212) + n^2*(2*A2112+2*A2211-A2222-A1111);
B3 = 2*r^2*A1212_p+6*r*A1212;
B4 = r^2*A1212;

%ode system
yp(1) = y(2);%yp(1)=f'(r)
yp(2) = y(3);%yp(2)=f''(r)
yp(3) = y(4);%yp(3)=f'''(r)
yp(4) = -(B0*y(1)+B1*y(2)+B2*y(3)+B3*y(4))/B4;%yp(4)=f''''(r)
yp(5) = W1*2*(-alpha^(-2) + alpha^2)/r;

%for a certain a, find b by r^2-a^2 = R^2-A^2 
%imposed by incompressible condition
%[R,Y] = ode45(odesystem, [a,b], ya1);
%[R,Y] = ode45(odesystem, [a,b], ya2);
%evaluate Y at b
%D(a) = bc1(yb1(1),yb1(2),yb1(3),yb1(4))*bc2(yb2(1),yb2(2),yb2(3),yb3(4))...
%      -bc1(yb2(1),yb1(2),yb1(3),yb1(4))*bc2(yb1(1),yb1(2),yb1(3),yb1(4));
