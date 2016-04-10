
A=1;
B=1.1;
C=1.5;
mu_i = 10;
mu_o = 1;
n=2;

g1_i=1.0;
g2_i=1.0;
g1_o=1.0;
g2_o=1.0;
epsilon=1e-6;

kstore=zeros(20,1);
Pstore=zeros(20,1);
for i=1:1:19
g1_i=sqrt(i*0.1);
g2_i=sqrt(i*0.1);
kl=0.3;
ku=1.5;
while kl<ku && kl+epsilon<ku
    km = (kl+ku)/2;
    a = km*g2_i*A;
    b=sqrt(g1_i*g2_i*(B^2-A^2)+a^2);
    c=sqrt(g1_o*g2_o*(C^2-B^2)+b^2);
    %calculate external pressure according to a and b
    %P = log(B/A) + 1/2*(A^2-a^2)*(1/b^2-1/a^2) + log(a/b);%needed to be updated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%first layer%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ya1 = zeros(5,1);
    ya1(1) = 0.1;
    ya1(2) = -0.1;
    %a^2*y0(3)+a*y0(2)+(n^2-1)*y0(1) = 0;
    ya1(3) = -(a*ya1(2)+(n^2-1)*ya1(1))/(a^2);
    ya1(4) = y4(ya1(1),ya1(2),ya1(3),A,a,mu_i,n,g1_i,g2_i,0); 
    ya1(5) = 0;
    [rr,yy] = ode45(@(r,y) yp(r,y,A,a,mu_i,n,g1_i,g2_i), [a,b], ya1); 
    yb1 = zeros(5,1);
    yb1(1) = yy(length(yy),1);
    yb1(2) = yy(length(yy),2);
    yb1(3) = yy(length(yy),3);
    yb1(4) = yy(length(yy),4);
    yb1(5) = yy(length(yy),5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ya2 = zeros(5,1);
    ya2(1) = 0.1;
    ya2(2) = 0.5;
    ya2(3) = -(a*ya2(2)+(n^2-1)*ya2(1))/(a^2);
    ya2(4) = y4(ya2(1),ya2(2),ya2(3),A,a,mu_i,n,g1_i,g2_i,0);
    ya2(5) = 0;
    [rr,yy] = ode45(@(r,y) yp(r,y,A,a,mu_i,n,g1_i,g2_i), [a,b], ya2); 
    yb2 = zeros(5,1);
    yb2(1) = yy(length(yy),1);
    yb2(2) = yy(length(yy),2);
    yb2(3) = yy(length(yy),3);
    yb2(4) = yy(length(yy),4);
    yb2(5) = yy(length(yy),5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%second layer%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %bc1(yc1(1),yc1(2),yc1(3),yc1(4),C,c,n,g1_o,g2_o)
    %bc2(yc2(1),yc2(2),yc2(3),yc2(4),C,c,n)
    ybb1 = zeros(5,1);
    ybb1(1) = yb1(1);
    ybb1(2) = yb1(2);
    Pnormal = bc1(yb1(1),yb1(2),yb1(3),yb1(4),B,b,mu_i,n,g1_i,g2_i)/b;%previously (g1_i*g2_i^2*B^2*b^3);
    Pshear  = mu_i*bc2(yb1(1),yb1(2),yb1(3),yb1(4),B,b,n);
    ybb1(3) = (Pshear/mu_o - (b*ybb1(2)+(n^2-1)*ybb1(1)))/b^2;
    %ybb1(4) = (Pnormal - ((2*g2^3*B^2*b^4+2*g1*g2^4*B^4*b^2)*ybb1(3)+(2*g2^3*B^2*b^3-2*g1*g2^4*n^2*B^4*b-g1*g2^4*B^4*b-g1*n^2*b^5)*ybb1(2)+...
    %          (n^2-1)*(2*g2^3*B^2*b^2-g1*g2^4*B^4)*ybb1(1))/(g1_o*g2_o^2*B^2*b^3) )/(g2^2*B^2);
    ybb1(4) = y4(ybb1(1),ybb1(2),ybb1(3),B,b,mu_o,n,g1_o,g2_o,Pnormal*b);     
    ybb1(5) = yb1(5);
    [rr,yy] = ode45(@(r,y) yp(r,y,B,b,mu_o,n,g1_o,g2_o), [b,c], ybb1);
    yc1 = zeros(4,1);
    yc1(1) = yy(length(yy),1);
    yc1(2) = yy(length(yy),2);
    yc1(3) = yy(length(yy),3);
    yc1(4) = yy(length(yy),4);
    P = -yy(length(yy),5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ybb2 = zeros(5,1);
    ybb2(1) = yb2(1);
    ybb2(2) = yb2(2);
    Pnormal = bc1(yb2(1),yb2(2),yb2(3),yb2(4),B,b,mu_i,n,g1_i,g2_i)/b;%previously(g1_i*g2_i^2*B^2*b^3);
    Pshear  = mu_i*bc2(yb2(1),yb2(2),yb2(3),yb2(4),B,b,n);
    ybb2(3) = (Pshear/mu_o - (b*ybb2(2)+(n^2-1)*ybb2(1)))/(b^2);
    ybb2(4) = y4(ybb2(1),ybb2(2),ybb2(3),B,b,mu_o,n,g1_o,g2_o,Pnormal*b);
    ybb2(5) = 0;
    [rr,yy] = ode45(@(r,y) yp(r,y,B,b,mu_o,n,g1_o,g2_o), [b,c], ybb2);
    yc2 = zeros(4,1);
    yc2(1) = yy(length(yy),1);
    yc2(2) = yy(length(yy),2);
    yc2(3) = yy(length(yy),3);
    yc2(4) = yy(length(yy),4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%evaluate Y at b
    %D = bc1(yc1(1),yc1(2),yc1(3),yc1(4),C,c,n,g1_o,g2_o)*bc2(yc2(1),yc2(2),yc2(3),yc2(4),C,c,n)...
    %      -bc1(yc2(1),yc2(2),yc2(3),yc2(4),C,c,n,g1_o,g2_o)*bc2(yc1(1),yc1(2),yc1(3),yc1(4),C,c,n);
    %D = - yc1(1)*bc2(yc2(1),yc2(2),yc2(3),yc2(4),C,c,n) + yc2(1)*bc2(yc1(1),yc1(2),yc1(3),yc1(4),C,c,n);
    D = -yc1(1)*yc2(2) + yc1(2)*yc2(1);
    if(D>0) kl=km;
    else if(D<0) ku=km;
        else break;
        end
    end
end
kstore(i)=ku;
Pstore(i)=P;
end
kstore
Pstore
