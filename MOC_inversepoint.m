function [up,vp,x2,y2] = MOC_inversepoint(x1,y1,u1,v1,x3,y3,u3,v3,x4,y4,theta4,To,R,gamma,del,Po)

%% initially we staet with pt. 2 and only the flow properties at pt. 2 are assumed to be
%  the same as that of pt. 3. 
u2 = u3 ; v2 = v3 ;
% now we start solving for lambda2 i.e lambda positive corrospond tpo L.R.C C+
%% first we will solve for lemda positive and lembda negative 
for i = 1:10  % usually take large number of iterations.
thetap1 = atan(v2/u2) ;% theta positive
thetap = (thetap1)*180/pi ;
modv2 =(u2.^2 + v2.^2).^0.5 ;    % modulus of v+
ao = (gamma * R * To)^0.5 ;
a2 = (ao.^2 - ((gamma - 1)/2)*modv2.^2).^0.5 ; % local acoustic speeed positive.
Mp = modv2./a2 ;
alphap1 = asin(1./Mp);
alphap = alphap1*180/pi ;
lambda2 = tan(thetap1+alphap1);
% solving for the linear equation 13 and putting oint 2 on that line.
m31 = (y1 - y3)/(x1 - x3) ;
% solving for location of point 2 ie. x2 and y2 
A = [1,-m31;1,-lambda2];
B = [y1 - m31*x1;y4 - lambda2*x4];
C  = linsolve(A,B);
y2 = C(1) ;
x2 = C(2) ;
% linear interpolation to solve for velocity at point 2 
f = (x2 - x1)/(x3 - x1) ;
u2i = u1 + f*(u3 - u1) ;
v2i = v1 + f*(v3 - v1) ;
u2 = u2i ;
v2 = v2i ;
end
% take the values of u2 and v2 that we obtained after above n iterations
% solve for all the positive flow parametres i.e C+ L.R.C
ICOR = 3 ;
V = zeros(2,ICOR);

for n = 1:ICOR
thetap1 = atan(v2/u2) ;% theta positive
thetap = (thetap1)*180/pi ;
modv2 =(u2.^2 + v2.^2).^0.5 ;    % modulus of v+
ao = (gamma * R * To)^0.5 ;
a2 = (ao.^2 - ((gamma - 1)/2)*modv2.^2).^0.5 ; % local acoustic speeed positive.
Mp = modv2./a2 ;
alphap1 = asin(1./Mp);
alphap = alphap1*180/pi ;
lambda2 = tan(thetap1+alphap1);

Q2 = (u2.^2 - a2.^2);
R2 = 2*u2*v2 - Q2*lambda2 ;
S2 = (del*(a2.^2)*v2)./y2 ;
T2 = S2*(x4-x2) + Q2*u2 + R2*v2 ;

%% solving for flow parameters at u4 and v4 i.e the required wll point fkow parametres.
A1 = [Q2,R2;tand(theta4),-1];
B1 = [T2;0];
C1 = linsolve(A1,B1) ;
u4p = C1(1) ; v4p = C1(2) ; % p stands for predictor/

u2 = (u2 + u4p)/2 ;
v2 = (v2 + v4p)/2;
s = [u2,v2]' ;
V(:,n) = V(:,n) + ([u4p v4p]') ; % predictor values.
end
V ;
up = V(1,ICOR) ;
vp = V(2,ICOR) ;

end