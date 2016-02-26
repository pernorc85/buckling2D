%for a certain a, find b by r^2-a^2 = R^2-A^2 
%imposed by incompressible condition
A=1;
B=1.2;
C=1.5;
mu_i = 5;%inner layer stiffness
mu_o = 1;%outer layer stiffness
k1=1;
k2=0.2;
g1_i = 1.0;%inner layer radical growth
g2_i = 1.0;%inner layer circumferential growth
g1_o = 1.0;%outer layer radical growth
g2_o = 1.0;%outer layer circumferential growth

epsilon=1e-6;
PI=3.141592654;
gamma=PI/18;

kstore=zeros(10,15);
Pstore=zeros(10,15);
for n=6:1:7
for i=31:1:40
g1_i=0.2*i;
%g2_i=0.2*i;
kl=0.3;
ku=0.7;
while kl<ku && kl+epsilon<ku
    km = (kl+ku)/2;
%for ii=1:1:9
    a = km*g2_i*A;
    b = sqrt( (B^2-A^2)*g1_i*g2_i + a^2);
    c = sqrt( (C^2-B^2)*g1_o*g2_o + b^2);
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
    ya1(4) = y4_aniso(ya1(1),ya1(2),ya1(3),A,a,n,mu_i,k1,k2,gamma,g1_i,g2_i);
    ya1(5) = 0;
    [rr,yy] = ode45(@(r,y) yp_aniso(r,y,A,a,mu_i,k1,k2,n,gamma,g1_i,g2_i), [a,b], ya1); 
    subplot(2,2,1);
    plot(rr,yy(:,1),rr,yy(:,2),rr,yy(:,3),rr,yy(:,4));
    legend('y1','y2','y3','y4');
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
    ya2(4) = y4_aniso(ya2(1),ya2(2),ya2(3),A,a,n,mu_i,k1,k2,gamma,g1_i,g2_i);
    ya2(5) = 0;
    [rr,yy] = ode45(@(r,y) yp_aniso(r,y,A,a,mu_i,k1,k2,n,gamma,g1_i,g2_i), [a,b], ya2); 
    subplot(2,2,3);
    plot(rr,yy(:,1));
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
    Pnormal = bc1_aniso(yb1(1),yb1(2),yb1(3),yb1(4),B,b,n,mu_i,k1,k2,gamma,g1_i,g2_i)/b;%previously it was /B^2/b^3
    Pshear  = mu_i*bc2(yb1(1),yb1(2),yb1(3),yb1(4),B,b,n);
    ybb1(3) = (Pshear/mu_o - (b*ybb1(2)+(n^2-1)*ybb1(1)))/b^2;
    ybb1(4) = Pnormal/g2_o^2/B^2 + y4(ybb1(1),ybb1(2),ybb1(3),B,b,n,g1_o,g2_o);
    ybb1(5) = yb1(5);
    [rr,yy] = ode45(@(r,y) yp(r,y,B,b,mu_o,n,g1_o,g2_o), [b,c], ybb1);
    subplot(2,2,2);
    plot(rr,yy(:,1),rr,yy(:,2),rr,yy(:,3),rr,yy(:,4));
    legend('y1','y2','y3','y4');
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
    Pnormal = bc1_aniso(yb2(1),yb2(2),yb2(3),yb2(4),B,b,n,mu_i,k1,k2,gamma,g1_i,g2_i)/b;%previously it was /B^2/b^2
    Pshear  = mu_i*bc2(yb2(1),yb2(2),yb2(3),yb2(4),B,b,n);
    ybb2(3) = (Pshear/mu_o - (b*ybb2(2)+(n^2-1)*ybb2(1)))/b^2;
    ybb2(4) = Pnormal/g2_o^2/B^2 + y4(ybb2(1),ybb2(2),ybb2(3),B,b,n,g1_o,g2_o);
    ybb2(5) = 0;
    [rr,yy] = ode45(@(r,y) yp(r,y,B,b,mu_o,n,g1_o,g2_o), [b,c], ybb2);
    subplot(2,2,4);
    plot(rr,yy(:,1),rr,yy(:,2),rr,yy(:,3),rr,yy(:,4));
    legend('y1','y2','y3','y4');
    yc2 = zeros(4,1);
    yc2(1) = yy(length(yy),1);
    yc2(2) = yy(length(yy),2);
    yc2(3) = yy(length(yy),3);
    yc2(4) = yy(length(yy),4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%evaluate Y at b
    D = bc1(yc1(1),yc1(2),yc1(3),yc1(4),C,c,n,g1_o,g2_o)*bc2(yc2(1),yc2(2),yc2(3),yc2(4),C,c,n)...
       -bc1(yc2(1),yc2(2),yc2(3),yc2(4),C,c,n,g1_o,g2_o)*bc2(yc1(1),yc1(2),yc1(3),yc1(4),C,c,n);
    %D = bc1(yc2(1),yc2(2),yc2(3),yc2(4),C,c,n,g1_o,g2_o)*yc1(1)...
    %     -bc1(yc1(1),yc1(2),yc1(3),yc1(4),C,c,n,g1_o,g2_o)*yc2(1);
    %D = - yc1(1)*yc2(2) + yc1(2)*yc2(1);
    if(D>0) kl=km;
    else if(D<0) ku=km;
        else break;
        end
    end 
end

kstore(i,n)=ku;
if(P>100)P=100;
elseif(P<-100)P=-100;
end
Pstore(i,n)=P;
end
end
kstore
Pstore
