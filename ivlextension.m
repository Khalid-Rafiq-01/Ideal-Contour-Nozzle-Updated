%{
  Code Description:
  Code is for the extension of I.V.L and is to be run after the file seuers method. It explains the
  pattern of flowfield emerging beyond the supersonic IVL. 
  Matrixes x,y,u,v each of size n1*n1(no. of points chosen on I.V.L) are generated
  Ist row of x,y,u,v contains flow properties on these points on I.V.L. Then an upper 
  triangular matrix is generated for each x,y,u,v that only requires function file 
  MOC_int_points. The generation of lower triangular
  matrix involves one function file but with two different sets of input
  such as one acts as an ordinary interior point solver while the other
  acts as an axis point solver. An if statement makes sure the respective function
  files are used for correctly solving the matrices. 

%}


x(1,1) = xivlbar(1);
y(1,1) = 10^(-9);    % choose a small value to avoid singularity.
u(1,1) = u(1);
v(1,1) = 0 ;

% Creating an x,y,u,v matrix with size n1*n1. First row is all about flow properties of I.V.L.

for i = 2:n1 
    x(1,i)=xivlbar(i) ;
    y(1,i)=y(i);
    u(1,i)=u(i);
    v(1,i)=0 ;
end

% Generation of the upper triangular matrix.

for j = 2:n1
    for k = j:n1
        [x(j,k),y(j,k),u(j,k),v(j,k)] = MOC_int_point(x(j-1,k-1),y(j-1,k-1),u(j-1,k-1),v(j-1,k-1),...
                                                      x(j-1,k),y(j-1,k),u(j-1,k),v(j-1,k),To,R,gamma,del,Po);
    end
end

% Generation of the lower triangular matrix.
for j = 2:n1
    for k = j:-1:2
        if j*(k-1)-j == 0
         [x(j,k-1),y(j,k-1),u(j,k-1),v(j,k-1)] =   MOC_symmetry(x(j,2),y(j,2),u(j,2),v(j,2),To,R,gamma,del,Po) ;
                                                  
             %{
               
               For a point on the axis of symmetry, interior point solver function
               can be used if one point is known and for the other point,
               the value of x and u are taken the same as the known point but the 
               values of y and v are taken as negatives of the known point.
                
             %}
        else
          [x(j,k-1),y(j,k-1),u(j,k-1),v(j,k-1)] = MOC_int_point(x(j-1,k-2),y(j-1,k-2),u(j-1,k-2),v(j-1,k-2),...
                                                      x(j,k),y(j,k),u(j,k),v(j,k),To,R,gamma,del,Po);
        end
    end
end

% Plotting the upper triangular matrix.
for j = 1:(n1-1) 
    for i = 1:(n1-(j-1)) 
    x1(i) = x(i,i+j-1);
    y1(i) = y(i,i+j-1) ;
    plot(x1,y1,'m','linewidth',2)
    end
end

hold on 
plot([x(1,n1),x(2,n1)],[y(1,n1),y(2,n1)],'m','linewidth',2)
hold on 

% Plotting the lower triangular matrix.

for j = 1:(n1-1) 
    for i = 1:(n1-(j-1)) 
    x2(i) = x(i+j-1,i);
    y2(i) = y(i+j-1,i) ;
    plot(x2,y2,'m','linewidth',2)
    end
end

hold on 
plot([x(n1,1),x(n1,2)],[y(n1,1),y(n1,2)],'m','linewidth',2)


% Solving for various Thermodynamic properties and Mach numbers.
 
 Modv = (u.^2 + v.^2).^0.5 ;
 ao = (gamma * R * To)^0.5 ;
 a = (ao.^2 - ((gamma - 1)/2)*Modv.^2).^0.5 ;
 Mivlx = Modv./a ;
 
 
 % Solving for thermodynamic properties.
 
 b = (gamma-1)/ 2 ;
 pivlx = ((1 + b*Mivlx.^2).^(-gamma/(2*b)))*Po ;
 tivlx = ((1 + b*Mivlx.^2).^(-1))*To ;
 rho = pivlx./(R*tivlx);
 

 

 

