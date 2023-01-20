%{
 
 This part is the continuation of Saeurs method and the IVL extension
 method.Here a circular arc is taken as the initial extension region and
 solution field is extended here by using the Method of Characteristics
 until a desired exit mach number is achieved at the axis of symmetry. Then
 the solution is checked by using the concept of consistency of mass flow
 rates and comparing the value of Area-Mach relationship given by the code and
 one provided by online NASA calclator. 

%}
%%
a = 2*length(x)-1 ;
rc = input('Enter the downstream radius of curvature: ') ;    % Radius of curvature of circular arc downstream of the throat.
b = n1 + 1 ;
theta = 0 ; % Initialization of theta parameter.
c = 0 ;
Md = input('Enter the exit mach number: ')  ;

while M<=Md   % Md is the desired Mach number.
 dtheta = 0.25; % The incremental value of theta that we want to creat the arc from.
 theta = theta + dtheta ;
  for i = b:b 
      for j = 1:a+c
          if (i*j)-i == 0 && i-n1-1==0 && a-j+c ~= 0 % First wall point on n1+1 line.
             [x(i,j),y(i,j)] = circle_arc(rc,theta,th) ;
             [u(i,j),v(i,j),xa,ya] = MOC_inversepoint(x(2,n1),y(2,n1),u(2,n1),v(2,n1),x(1,n1),...
                                                    y(1,n1),u(1,n1),v(1,n1),...
                                                         x(i,j),y(i,j),theta,To,R,gamma,del,Po);
                                                    
                       if xa>x(1,n1) && xa<x(2,n1) && ya<y(1,n1) && ya>y(2,n1)  % If Ist line does not fail.
                          [x(i,j+1),y(i,j+1),u(i,j+1),v(i,j+1)] = MOC_int_point(x(2,n1),y(2,n1),u(2,n1),...
                                                                        v(2,n1), x(i,j),y(i,j),u(i,j),...
                                                                            v(i,j),To,R,gamma,del,Po);
                             
                                hold on
                                
                                
                           plot([x(i,j) x(i,j+1)],[y(i,j) y(i,j+1)],'linewidth',2);
                              hold on 
                           plot([x(1,n1) x(i,j)],[y(1,n1) y(i,j)],'linewidth',2);
                              hold on 
                           plot([x(2,n1) x(i,j+1)],[y(2,n1) y(i,j+1)],'linewidth',2);
                           
                       elseif xa<x(1,n1) || xa>x(2,n1) || ya>y(1,n1) || ya<y(2,n1) % If Ist line fails.
                            [u(i,j),v(i,j),xa,ya] = MOC_inversepoint(x(3,n1),y(3,n1),u(3,n1),v(3,n1),...
                                                                        x(2,n1),y(2,n1),u(2,n1),v(2,n1),...
                                                                             x(i,j),y(i,j),theta,To,R,gamma,del,Po);
                                                                         
                            [x(i,j+1),y(i,j+1),u(i,j+1),v(i,j+1)] = MOC_int_point(x(3,n1),y(3,n1),u(3,n1),...
                                                                               v(3,n1), x(i,j),y(i,j),u(i,j),...
                                                                                          v(i,j),To,R,gamma,del,Po);
                                                                                      
                               hold on
                               
                           plot([x(i,j) x(i,j+1)],[y(i,j) y(i,j+1)],'linewidth',2);
                               hold on 
                           plot([x(1,n1) x(i,j)],[y(1,n1) y(i,j)],'linewidth',2);
                              hold on 
                           plot([x(3,n1) x(i,j+1)],[y(3,n1) y(i,j+1)],'linewidth',2);
                               c = c-1 ;
                       end
                       
                           
          elseif  (i*j)-i == 0  && a-j+c ~= 0 % General wall point.
             [x(i,j),y(i,j)] = circle_arc(rc,theta,th) ;
             [u(i,j),v(i,j),xb,yb] = MOC_inversepoint(x(i-1,2),y(i-1,2),u(i-1,2),v(i-1,2),x(i-1,1),...
                                                                           y(i-1,1),u(i-1,1),v(i-1,1),...
                                                                             x(i,j),y(i,j),theta,To,R,gamma,del,Po); 
                        if xb>x(i-1,1) && xb<x(i-1,2) && yb<y(i-1,1) && yb>y(i-1,2) % no fail condition.
                                  [x(i,j+1),y(i,j+1),u(i,j+1),v(i,j+1)] = MOC_int_point(x(i-1,2),y(i-1,2),...
                                                                                           u(i-1,2),v(i-1,2),x(i,j),...
                                                                                                 y(i,j),u(i,j),v(i,j),...
                                                                                                      To,R,gamma,del,Po) ;
                                       hold on 
                                   plot([x(i,j) x(i,j+1)],[y(i,j) y(i,j+1)],'linewidth',2)     % Plots Right running char.
                                       hold on 
                                   plot([x(i-1,2) x(i,j+1)],[y(i-1,2) y(i,j+1)],'linewidth',2) % Plots bottom Left running char.
                                       hold on 
                                   plot([x(i-1,1) x(i,j)],[y(i-1,1) y(i,j)],'linewidth',2)     % Plots wall points LRC
                                       
                        elseif xb<x(i-1,1) || xb>x(i-1,2) || yb>y(i-1,1) || yb<y(i-1,2)        % Fail condition and wall point.
                                       [u(i,j),v(i,j),xb,yb] = MOC_inversepoint(x(i-1,3),y(i-1,3),u(i-1,3),...
                                                                                   v(i-1,3),x(i-1,2),y(i-1,2),...
                                                                                      u(i-1,2),v(i-1,2),x(i,j),y(i,j),...
                                                                                                    theta,To,R,gamma,del,Po);
                                                                                                
                                       [x(i,j+1),y(i,j+1),u(i,j+1),v(i,j+1)] = MOC_int_point(x(i-1,3),y(i-1,3),...
                                                                                               u(i-1,3),v(i-1,3),x(i,j),...
                                                                                                 y(i,j),u(i,j),v(i,j),To,...
                                                                                                           R,gamma,del,Po);
                                         hold on 
                                        plot([x(i,j) x(i,j+1)],[y(i,j) y(i,j+1)],'b','linewidth',2) % plots fail rrc
                                         hold on 
                                        plot([x(i-1,3) x(i,j+1)],[y(i-1,3) y(i,j+1)],'linewidth',2) % plots bottom lrc
                                         hold on 
                                        plot([x(i-1,1) x(i,j)],[y(i-1,1) y(i,j)],'linewidth',2) % plots wall points lrc
                                                                                 
                                                                       c = c-1;        
                        end
        % Not wall point ,Ist line,upper part and no failing of Ist line.                       
          elseif (i*j)-i ~= 0 && i-n1-1 == 0  && i-j-1 > 0 && a-j+c ~= 0 &&...
                                        xa>x(1,n1) && xa<x(2,n1) && ya<y(1,n1) && ya>y(2,n1) 
              
                           [x(i,j+1),y(i,j+1),u(i,j+1),v(i,j+1)] = MOC_int_point(x(j+1,n1),y(j+1,n1),...
                                                                               u(j+1,n1),v(j+1,n1), x(i,j),y(i,j),...
                                                                                       u(i,j),v(i,j),To,R,gamma,del,Po);
                                                   hold on
                                         plot([x(i,j) x(i,j+1)],[y(i,j) y(i,j+1)],'linewidth',2)
                                                   hold on 
                                         plot([x(j+1,n1) x(i,j+1)],[y(j+1,n1) y(i,j+1)],'linewidth',2)
                                                             
       %  Not wall point ,Ist line,lower part and Ist line does not fail. 
          elseif  (i*j)-i ~= 0 && i-n1-1 == 0  && i-j-1 <= 0 && a-j+c ~= 0 &&...
                                                    xa>x(1,n1) && xa<x(2,n1) && ya<y(1,n1) && ya>y(2,n1)
                           [x(i,j+1),y(i,j+1),u(i,j+1),v(i,j+1)] = MOC_int_point(x(i-1,2*n1-1-j),y(i-1,2*n1-1-j),...
                                                                                    u(i-1,2*n1-1-j),v(i-1,2*n1-1-j),...
                                                                                          x(i,j),y(i,j),u(i,j),v(i,j),...
                                                                                                      To,R,gamma,del,Po);
                                                   hold on 
                                        plot([x(i,j) x(i,j+1)],[y(i,j) y(i,j+1)],'linewidth',2)
                                                   hold on 
                                        plot([x(i-1,2*n1-1-j) x(i,j+1)],[y(i-1,2*n1-1-j) y(i,j+1)],'linewidth',2)
               
                % Not wall point ,Ist line,upper part and failing of Ist line.                       
          elseif (i*j)-i ~= 0 && i-n1-1 == 0  && i-j-2 > 0 && a-j+c ~= 0 &&...
                                          (xa<x(1,n1) || xa>x(2,n1) || ya>y(1,n1) || ya<y(2,n1)) 
              
                            [x(i,j+1),y(i,j+1),u(i,j+1),v(i,j+1)] = MOC_int_point(x(j+2,n1),y(j+2,n1),...
                                                                                    u(j+2,n1),v(j+2,n1), x(i,j),...
                                                                                       y(i,j),u(i,j),v(i,j),To,R,...
                                                                                                         gamma,del,Po);
                                                   hold on
                                        plot([x(i,j) x(i,j+1)],[y(i,j) y(i,j+1)],'linewidth',2)
                                                   hold on 
                                        plot([x(j+2,n1) x(i,j+1)],[y(j+2,n1) y(i,j+1)],'linewidth',2)
                
             %  Not wall point ,Ist line,lower part and Ist line fails.                                                  
          elseif  (i*j)-i ~= 0 && i-n1-1 == 0  && i-j-2 <= 0 && a-j+c ~= 0 &&...
                                           (xa<x(1,n1) || xa>x(2,n1) || ya>y(1,n1) || ya<y(2,n1)) 
               
                            [x(i,j+1),y(i,j+1),u(i,j+1),v(i,j+1)] = MOC_int_point(x(i-1,2*n1-2-j),...
                                                                      y(i-1,2*n1-2-j),u(i-1,2*n1-2-j),...
                                                                             v(i-1,2*n1-2-j),x(i,j),y(i,j),...
                                                                                     u(i,j),v(i,j),To,R,gamma,del,Po);
                                                   hold on 
                                       plot([x(i,j) x(i,j+1)],[y(i,j) y(i,j+1)],'linewidth',2)
                                                   hold on 
                                       plot([x(i-1,2*n1-2-j) x(i,j+1)],[y(i-1,2*n1-2-j) y(i,j+1)],'linewidth',2)
               
             % Ist Line symmetry point and the Ist line fails.
          elseif  (i*j)-i ~= 0 && i-n1-1 == 0  && i-j-2 <= 0 && a-j+c == 0 &&...
                                                          (xa<x(1,n1) || xa>x(2,n1) || ya>y(1,n1) || ya<y(2,n1))
                                  [x(i,j+1),y(i,j+1),u(i,j+1),v(i,j+1)] =  MOC_symmetry( x(i,j),y(i,j),u(i,j),...
                                                                                          v(i,j),To,R,gamma,del,Po);
                                                   hold on 
                                       plot([x(i,j) x(i,j+1)],[y(i,j) y(i,j+1)],'linewidth',2)
          break
          
             % Not a wall point, not first line, and no fail.                                                 
          elseif (i*j)-i ~= 0 && i-n1-1 ~= 0 && xb>x(i-1,1) && xb<x(i-1,2) &&...
                                                          yb<y(i-1,1) && yb>y(i-1,2) && a-j+c ~= 0
                                  [x(i,j+1),y(i,j+1),u(i,j+1),v(i,j+1)] = MOC_int_point(x(i-1,j+1),y(i-1,j+1),...
                                                                                            u(i-1,j+1),v(i-1,j+1),...
                                                                                             x(i,j),y(i,j),u(i,j),v(i,j),...
                                                                                                         To,R,gamma,del,Po);
                                                   hold on
                                      plot([x(i,j) x(i,j+1)],[y(i,j) y(i,j+1)],'linewidth',2)
                                                   hold on 
                                      plot([x(i-1,j+1) x(i,j+1)],[y(i-1,j+1) y(i,j+1)],'linewidth',2)
                                                                
              % Not a wall point, not first line, and fails.                                                     
          elseif (i*j)-i ~= 0 && i-n1-1 ~= 0 && a-j+c ~= 0 && (xb<x(i-1,1) ||...
                                                         xb>x(i-1,2) || yb>y(i-1,1) || yb<y(i-1,2))
                                  [x(i,j+1),y(i,j+1),u(i,j+1),v(i,j+1)] = MOC_int_point(x(i-1,j+2),y(i-1,j+2),...
                                                                                 u(i-1,j+2),v(i-1,j+2),x(i,j),y(i,j),...
                                                                                          u(i,j),v(i,j),To,R,gamma,del,Po);
                                                   hold on 
                                      plot([x(i,j) x(i,j+1)],[y(i,j) y(i,j+1)],'b','linewidth',2)
                                                   hold on 
                                      plot([x(i-1,j+2) x(i,j+1)],[y(i-1,j+2) y(i,j+1)],'linewidth',2)
                 
          elseif a-j+c == 0 % approached line of symmetry cdtn.
                                  [x(i,j+1),y(i,j+1),u(i,j+1),v(i,j+1)] =  MOC_symmetry( x(i,j),y(i,j),u(i,j),...
                                                                                         v(i,j),To,R,gamma,del,Po);
                                                   hold on 
                                     plot([x(i,j) x(i,j+1)],[y(i,j) y(i,j+1)],'linewidth',2)
             
             break
              
          end
      end
      
           c = c+1 ;
           b = b+1 ;
         
  end
   modv =(u.^2 + v.^2).^0.5 ;                    % Magnitude of speed of air at different nodes.
   a1 = (ao.^2 - ((gamma - 1)/2)*modv.^2).^0.5 ; % Relative acoustic speed.
   M = modv./a1 ;
   
end   

 %% Calculation of area ratio using concept of mass flow rates.
 
  [a,b] = size(x);
  f = (gamma-1)/ 2 ; 
  mi = M(end,nnz(M(end,:))) 
  pl = ((1 + f*mi.^2).^(-gamma/(2*f)))*Po ;
  tl = ((1 + f*mi.^2).^(-1))*To ;
  rhol = pl./(R*tl);   % Density of the final solution point on the axis of symmetry.
  ui = u(end,nnz(u(end,:)));
  vi = v(end,end);
  mvi = modv(end,end); % Mod of velocity of sol. point on the axis of symmetry. 
  aratio = (Mdotnumeric/(pi*rhol*ui))/((th)^2) % Derived from conservation of mdot.(Refer to G.V.R Rao's paper on optimal contour)
  
  %% Matrix solution of P,T,Rho, and modulus of speed.
  Pm = ((1 + f*M.^2).^(-gamma/(2*f)))*Po ;
  Tm = ((1 + f*M.^2).^(-1))*To ;  
  Rho = Pm./(R*Tm);
  
%% Plotting for the end point.
alpha = asind(1/mi) ;
xl = x(end,nnz(x(end,:))) ;
yl = y(end,nnz(y(end,:))) ;
h = (aratio)^0.5 * th;
modLf = h*mi ;
xf = modLf*cosd(alpha) + xl ; % Represents the lenght of the nozzle.
yf = h ;
hold on 
plot([xl,xf],[yl,yf],'-k','linewidth',2)
hold on 
plot(xf,yf,'*','markersize',3) ;

  