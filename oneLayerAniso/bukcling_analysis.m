%for a certain a, find b by r^2-a^2 = R^2-A^2 
%imposed by incompressible condition
%A=0.3;
B=1;
%n=8;
PI=3.141592654;
gamma=PI/180*60;

kstore=zeros(20,15);
Pstore=zeros(20,15);
%Dstore=zeros(20,20);
for n=2:1:7
for i=1:1:5
    
A = i*0.05*B;
kl = 0.3;
ku = 0.99;
epsilon = 1e-6;
while kl<ku && kl+epsilon<ku
    km = (kl+ku)/2;
    a = km*A;
    b=sqrt(B^2-A^2+a^2);
    %calculate external pressure according to a and b
    P = log(B/A) + 1/2*(A^2-a^2)*(1/b^2-1/a^2) + log(a/b);%confirmed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ya1 = zeros(4,1);
    ya1(1) = 0.1;
    ya1(2) = -0.1;
    %A^4*a^3*y0(4)+(2*A^2*a^4+2*A^4*a^2)*y0(3)+(2*A^2*a^3-2*n^2*A^4*a-A^4*a-n^2*a^5)*y0(2)+(n^2-1)*(2*A^2*a^2-A^4)*y0(1) = 0;
    %a^2*y0(3)+a*y0(2)+(n^2-1)*y0(1) = 0;
    ya1(3) = -(a*ya1(2)+(n^2-1)*ya1(1))/(a^2);
    ya1(4) = y4(ya1(1),ya1(2),ya1(3),A,a,n,gamma);%-((2*A^2*a^4+2*A^4*a^2)*ya1(3)+(2*A^2*a^3-2*n^2*A^4*a-A^4*a-n^2*a^5)*ya1(2)+(n^2-1)*(2*A^2*a^2-A^4)*ya1(1))/(A^4*a^3);
    [rr,yy] = ode45(@(r,y) yp(r,y,A,a,n,gamma), [a,b], ya1);
    subplot(2,1,1);
    plot(rr,yy(:,1),rr,yy(:,2),rr,yy(:,3),rr,yy(:,4));
    legend('y1','y2','y3','y4');
    yb1 = zeros(4,1);
    yb1(1) = yy(length(yy),1);
    yb1(2) = yy(length(yy),2);
    yb1(3) = yy(length(yy),3);
    yb1(4) = yy(length(yy),4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ya2 = zeros(4,1);
    ya2(1) = 0.1;
    ya2(2) = 0.2;
    ya2(3) = -(a*ya2(2)+(n^2-1)*ya2(1))/(a^2);
    ya2(4) = y4(ya2(1),ya2(2),ya2(3),A,a,n,gamma);%-((2*A^2*a^4+2*A^4*a^2)*ya2(3)+(2*A^2*a^3-2*n^2*A^4*a-A^4*a-n^2*a^5)*ya2(2)+(n^2-1)*(2*A^2*a^2-A^4)*ya2(1))/(A^4*a^3);
    [rr,yy] = ode45(@(r,y) yp(r,y,A,a,n,gamma), [a,b], ya2);  
    subplot(2,1,2);
    plot(rr,yy(:,1));
    yb2 = zeros(4,1);
    yb2(1) = yy(length(yy),1);
    yb2(2) = yy(length(yy),2);
    yb2(3) = yy(length(yy),3);
    yb2(4) = yy(length(yy),4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%evaluate Y at b
%D = yb1(1)*bc2(yb2(1),yb2(2),yb2(3),yb2(4),B,b,n)...
%          -yb2(1)*bc2(yb1(1),yb1(2),yb1(3),yb1(4),B,b,n);
D = bc1(yb1(1),yb1(2),yb1(3),yb1(4),B,b,n,gamma)*bc2(yb2(1),yb2(2),yb2(3),yb2(4),B,b,n)...
    -bc1(yb2(1),yb2(2),yb2(3),yb2(4),B,b,n,gamma)*bc2(yb1(1),yb1(2),yb1(3),yb1(4),B,b,n);
    if(D>0) kl=km;
    else if(D<0) ku=km;
        else break;
        end
    end     
end

kstore(i,n)=ku;
Pstore(i,n)=P;
end
end
kstore
Pstore
