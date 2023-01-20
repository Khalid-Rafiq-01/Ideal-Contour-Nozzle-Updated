%{

  This code is used to find the flow properties of a node lying on the axis
  of symmetry. Internal point solver function can also be used to solve the 
  flow properties this node as the flow properties of the node above it are known
  and the 'x' and 'u' of the other point can be taken same while 'y' and
  'v' are taken as opposites of this.

%}

function [xp,yp,up,vp] = MOC_symmetry(x1,y1,u1,v1,To,R,gamma,del,Po)

ICOR = 3 ;
S = zeros(4,ICOR);
yp = 10^(-9) ;
vp = 0 ;
for n = 1: ICOR 
%% In this case a L.R.C will be used to join point 1 with point 4 which lies on the line of symmetry.
 % calculating the value of lambda1
 
thetan1 = atan(v1/u1) ;          % Theta positive
thetan = (thetan1)*180/pi ;
modv1 =(u1.^2 + v1.^2).^0.5 ;    % Modulus of v+
ao = (gamma * R * To)^0.5 ;
a1 = (ao.^2 - ((gamma - 1)/2)*modv1.^2).^0.5 ; % Local acoustic speeed positive.
Mn = modv1./a1 ;
alphan1 = asin(1./Mn);
alphan = alphan1*180/pi ;
lambda1 = tan(thetan1-alphan1);

xp = x1 - y1/lambda1 ;            % x coordinate of the point 4 lying on axis of symmetry.

%% now we will solve for the u component of velocity at point 4

Q1 = (u1.^2 - a1.^2);
R1 = 2*u1*v1 - Q1*lambda1 ;
S1 = (del*(a1.^2)*v1)./y1 ;
T1 = S1*(xp-x1) + Q1*u1 + R1*v1 ;

up = T1/Q1 ;                      % velocity at point 4


%% for the corrector, taking the average value of the properties 
u1 = (u1 + up)/2 ;
v1 = (v1 + vp)/2 ;
 S(:,n) = S(:,n) + ([xp yp up vp]') ; % Predictor values. The first element is x location and second is u 
 end
 S ;
xp = S(1,ICOR) ;
yp = S(2,ICOR) ;
up = S(3,ICOR) ;
vp = S(4,ICOR) ;
end
