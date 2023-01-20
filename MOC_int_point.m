%{ 
 
 This function is used to solve x,y,u,v of an unknown point 
 which is formed by the intersection of one point (x1,y1,u1,v1) and
 another point (x2,y2,u2,v2). Please make sure that through point
 1(x1,y1,u1,v1) a Left running chracteristic passes and through point 2(x2,y2,u2,v2)
 a Right running characteristic passes.Use of modified Euler corrector predictor
 algorithm has been made for generating and solving set of linear finite difference
 equations. only two correctors and used(ICOR = 3). For a corrector level of N
 use ICOR = (N+1). Here 3 stages of corrector are used.
(xp,yp,up,vp) are the solution values.

%}

function [xp,yp,up,vp] = MOC_int_point(x1,y1,u1,v1,x2,y2,u2,v2,To,R,gamma,del,Po)

ICOR = 3 ;        % Stages of corrector = 2
P = zeros(4,ICOR);
  for n = 1:ICOR

%% First we will solve for lemda positive and lembda negative 
thetan1 = atan(v1/u1) ;                          % Theta negative
thetan = (thetan1)*180/pi ;
modv1 =(u1.^2 + v1.^2).^0.5 ;                    % Modulus of v-
ao = (gamma * R * To)^0.5 ;
a1 = (ao.^2 - ((gamma - 1)/2)*modv1.^2).^0.5 ;   % Local acoustic speeed positive.
Mn = modv1/a1 ;
if Mn<=1
    xp = x1 ;
    thetap1 = atan(v2/u2) ;          % Theta positive
    thetap = (thetap1)*180/pi ;
    modv2 =(u2.^2 + v2.^2).^0.5 ;    % Modulus of v+
    ao = (gamma * R * To)^0.5 ;
    a2 = (ao.^2 - ((gamma - 1)/2)*modv2.^2).^0.5 ; % Local acoustic speeed positive.
    Mp = modv2./a2 ;
    alphap1 = asin(1./Mp);
    alphap = alphap1*180/pi ;
    lambda2 = tan(thetap1-alphap1);
    yp = y2 + lambda2*(x1-x2);
    
    %% Solving for the flow characteristics ;
    
    Q2 = (u2.^2 - a2.^2);
    R2 = 2*u2*v2 - Q2*lambda2 ;
    S2 = (del*(a2.^2)*v2)./y2 ;
    T2 = S2*(xp-x2) + Q2*u2 + R2*v2 ;

  
    up = u2;
    vp = (T2 - Q2*(up))/R2 ;
    
elseif Mn>1
    
alphan1 = asin(1./Mn);
alphan = alphan1*180/pi ;
lambda1 = tan(thetan1+alphan1);

%% Now we will solve for lembda positive.

thetap1 = atan(v2/u2) ;          % Theta positive
thetap = (thetap1)*180/pi ;
modv2 =(u2.^2 + v2.^2).^0.5 ;    % Modulus of v+
ao = (gamma * R * To)^0.5 ;
a2 = (ao.^2 - ((gamma - 1)/2)*modv2.^2).^0.5 ; % Local acoustic speeed positive.
Mp = modv2./a2 ;
alphap1 = asin(1./Mp);
alphap = alphap1*180/pi ;
lambda2 = tan(thetap1-alphap1);
 
%% Solve for coordinates of the solution point :
A = [1,-lambda2;1,-lambda1];
B = [y2 - lambda2*x2 ; y1 - lambda1*x1];
L = linsolve(A,B) ; % Location of the solution point [y,x]'
xp = L(2) ;   yp = L(1) ;

%% Solving for the flow characteristics ;

Q2 = (u2.^2 - a2.^2);
R2 = 2*u2*v2 - Q2*lambda2 ;
S2 = (del*(a2.^2)*v2)./y2 ;
T2 = S2*(L(2)-x2) + Q2*u2 + R2*v2 ;

Q1 = (u1.^2 - a1.^2);
R1 = 2*u1*v1 - Q1*lambda1 ;
S1 = (del*(a1.^2)*v1)./y1 ;
T1 = S1*(L(2)-x1) + Q1*u1 + R1*v1 ;

C = [Q2,R2;Q1,R1];
D = [T2;T1];

F = linsolve(C,D); % First the u comp. then the v comp. 
up = F(1) ;   vp = F(2) ;
end

% Algorithm for the corrector 
% Taking the average of the properties at the initial point and the predictor point. 
 
x1 = (x1 + xp)/2 ;
y1 = (y1 + yp)/2 ;
u1 = (u1 + up)/2  ;
v1 = (v1 + vp)/2 ;

x2 = (x2 + xp)/2 ;
y2 = (y2 + yp)/2 ;
u2 = (u2 + up)/2  ;
v2 = (v2 + vp)/2 ;

    P(:,n) = P(:,n) + ([xp yp up vp]') ; % Predictor values.
   
 end
P ; % Displays the matrix in column is the predictor, the rest are correctors and iteration of correctors.
xp = P(1,ICOR);
yp = P(2,ICOR);
up = P(3,ICOR);
vp = P(4,ICOR);
end








