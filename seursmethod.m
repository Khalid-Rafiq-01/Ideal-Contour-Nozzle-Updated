%{ 
  Code description:  
    This piece of code will generate the I.V.L(Initial value line) and the sonic line in the nozzle throat.
    Seuers Method is used to solve the transonic equation in nozzle throat.
    Use a throat height (th) = 1, to keep it dimensionless. 
    Required output values on I.V.L are xivlbar, y , u (Please note on I.V.L, value of v = 0)
    For a throat height of unity, the minimum value of radius of curvature downstream
    is 2 (ROCth = 2), as seuers method works for Normalized throat ROC > 2.
    
%}

clear all
format long

%% Gas properties. 
R = input('Enter the value if gas constant(R): ');         %Default = 287
gamma = input('Enter the value of gamma: ');               %Default = 1.4

%% Stagnation properties.
To = input('Enter the stagnation temperature:  ');      %Default = 800
Po = input('Enter the stagnation pressure:  ');         %Default = 7*10^6

%% Geometric properties
del = 1;%input('Enter the value of delta : '); % del = 0, for planar case. and 1 for axis-symmetric. We have axis symmetric case. 
th = input('Enter the height of throat: ' );                % throat height. Default = 0.025
ROCth = input('Enter the radius of curvature of throat: '); % Default = 0.050
n1 =  input('Enter the number of points on the Sonic and I.V.L,choose an odd number: ');
dl = th/(n1-1);                                             % height of each segment. 

%% sonic and initial value line calculations.
eta = -(th/(2*(3+del)))*((th*(gamma + 1)*(1+del))/ROCth)^0.5 ; % sonic line cuts the LOS.
alpha = ((1+del)/((gamma+1)*ROCth*th))^0.5 ;                   % coefficient of axial velocity pertubration.
xs = @(y) - ((gamma+1)*alpha*y.^2)/(2*(del+1)) ;               % the sonic line .
xivl = @(y) - ((gamma+1)*alpha*y.^2)/(2*(del+3));              % the initial value line.
b = (gamma-1)/ 2 ;
y = 0:dl:th ;                                                  % diviting throat height into segments,each of height dl.
Xs = xs(y);                                                    % conversion to mm.
Xivl = xivl(y);
xsbar = Xs - eta ;                                             % values of x coordinates of n1 points on the sonic line.
xivlbar = Xivl - eta ;                                         % values of x coordiates of n1 points on I.V.L 

%%
%{ 
 For plotting the Sonic line.
 
 h1 = plot(xsbar,y,'b','linewidth',2);
 h3 = get(h1,'Parent');
 set(h3,'Fontsize',13,'Linewidth',2);
 xlim([-0.40 xivlbar(1)+0.6]);
 ylim([0 th+0.15]);
 title('Sonic line in nozzle throat');
 xlabel('Axial location of the nozzle')
 ylabel('Radial location of the Nozzle')
 hold on
 plot(xsbar,y,'v');
%}

%%

   %For plotting the I.V.L 
figure(1)
   h2 = plot(xivlbar,y,'linewidth',2);
   h4 = get(h2,'Parent');
   set(h4,'Fontsize',13,'Linewidth',2);
   plot(xivlbar,y)
   xlim([0 xivlbar(1)+ 15.5]);
   ylim([0 th+ 5.25]);
   title('Extension of Initial value line in nozzle throat');
   xlabel('Axial location of the nozzle')
   ylabel('Radial location of the Nozzle')
   hold on 
%  legend('Sonic line','Initial value line') % If printing both the lines together, uncomment/

%}

%% Properties along the initial value line.
astar = ((2*gamma*R*To)/(gamma+1))^0.5;
udash = alpha*Xivl+((gamma+1)*(alpha^2)*(y.^2))/(2*(1+del)) ; % non dimensional pertubration x field
u = astar*(1+udash);
ao = (gamma * R * To)^0.5 ;
a = (ao.^2 - ((gamma - 1)/2)*u.^2).^0.5 ;
v = zeros(1,1*n1);
M = u./a ;
p = ((1 + b*M.^2).^(-gamma/(2*b)))*Po ;
t = ((1 + b*M.^2).^(-1))*To ;
rho = p./(R*t);

%{
  % This block shows the variation of Mach number on the I.V.L along the
  height of the throat.

  h5 = plot(M,y,'*','markersize',2);
  h6 = get(h5,'Parent');
  set(h6,'Fontsize',13,'Linewidth',2);
  xlim([0.8 M(n1)+0.5]);
  ylim([0 th+0.05]);
  title('Variation of Mach number along the radial direction along I.V.L');
  xlabel('Mach number')
  ylabel('Height of nozzle')
 
%}
 
%%
%{
  %variation of static pressure number along the height of throat along the I.V.L 
 
  h7 = plot(p,y,'-o','markersize',2);
  h8 = get(h7,'Parent');
  set(h8,'Fontsize',13,'Linewidth',2);
  xlim([2.5*10^6 p(1)+500000]);
  ylim([0 th+0.05]);
  title('Variation of pressure along the radial direction along I.V.L');
  xlabel('Static pressure')
  ylabel('Height of nozzle')

%}


%% calculating the mass flow rate (Simpsons method is used).

if del == 1  % axis-symmetric nozzle.
     C1 = (2*pi*dl)/3 ;
     S = rho.*u.*y;
     st4 = 0 ;
     st2 = 0;
         for i = 2:2:n1
           st4 = st4 + 4*S(i);
         end
               for i = 3:2:n1-1
                    st2 = st2 +2*S(i);
               end
   Mdotnumeric = C1*(S(1)+st4+st2+S(n1))
   Mdotanalytic = rho(1)*astar*pi*th^2 ;
   Cd = Mdotnumeric/Mdotanalytic
   
elseif del==0  % Condition for planar nozzle.
    w = input('Enter the width of planar nozzle:');
    C1 = (w*dl)/3 ;
    S = rho.*u';
    st4 = 0 ;
     st2 = 0;
         for i = 2:2:n1
           st4 = st4 + 4*S(i);
         end
               for i = 3:2:n1-1
                    st2 = st2 +2*S(i);
               end
   Mdotnumeric = C1*(S(1)+st4+st2+S(n1)) % only for half of the nozzle.
   Mdotanalytic = rho(1)*astar*w*th*2 ;
    Cd = (2*Mdotnumeric)/Mdotanalytic    % Coefficient of discharge.
end
