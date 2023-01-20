%{
  
  Here we define the final section of the problem. i.e The turning section.
  This section defines the contour responsible for turning the fluid in the
  axial direction and here it is developed by keeping the consistency of
  mass flow rate into consideration.

%}
n = input('Enter the number of parts to divide the final solution line into : ') ;
Le = modLf/n ;
[p,q] = size(x) ;
k = zeros(p,q) ; % Basic function term to calculate the mass flow rate. f(y)

%% Initialising the matrix.

xr = zeros(n,q) ;
yr =   zeros(n,q) ;
ur =   zeros(n,q) ;
vr =   zeros(n,q) ;
mv =   zeros(n,q) ;
mr =   zeros(n,q) ;
rhor = zeros(n,q) ;
kr =   zeros(n,q) ;
alr =  zeros(n,q) ; % Local acoustic speed.
pm =   zeros(n,q) ;
tm =   zeros(n,q) ;

%% Generating the matrix of extension red line.
for i = 1:n
    xr(i) = Le*i*cosd(alpha) + xl ;
    yr(i) = Le*i*sind(alpha) ;
    ur(i) = ui ;
    vr(i) = vi ;
    mv(i) = mvi ;
    mr(i) = mi ;
    rhor(i) = rhol ;
    alr(i) = a1(end,end) ; 
    pm(i) = pl ;
    tm(i) = tl ;
    
      if i == 1
          kr(i) = (rhor(i)*mv(i)*sind(atand((yr(i)/(xr(i)-xl)) - ...
                                               atand(vr(i)/ur(i)))) *yr(i))/...
                                                     (sind(atand(yr(i)/(xr(i)-xl)))) ;
          mline(i) = kr(i,1)*pi*yr(i,1) ;
          
      elseif i >1
          
          kr(i) = (rhor(i)*mv(i)*sind(atand(((yr(i)-yr(i-1))/(xr(i)-xr(i-1)))...
                                  - atand(vr(i)/ur(i))))*yr(i))/(sind(atand(((yr(i)-yr(i-1))...
                                                                            /(xr(i)-xr(i-1)))))) ;
          mline(i) = (kr(i,1)+kr(i-1,1))*pi*(yr(i,1) - yr(i-1,1)) ;
      end
   
    
end

%% concatinating the solution ine matrix and the original matrixes ...

x =    [x;xr] ;
y =    [y;yr] ;
u =    [u;ur] ;
v =    [v;vr] ;
Rho =  [Rho;rhor] ;
M =    [M;mr] ;
modv = [modv;mv] ;
k =    [k;kr] ;
al =   [a1;alr] ;
Pm =   [Pm;pm] ;
Tm =   [Tm;tm] ;


%%
%{
    Computing for interior points and also the mass flow rates.
    The first turning line follows a differnet nodal arrangement than the
    rest. That is accomidated into the code by using the if statement!
 
%}
 for i = p:p+n
     for j = 2:q 
         if i == p 
             [x(i+1,j),y(i+1,j),u(i+1,j),v(i+1,j)] = MOC_int_point(x(i,end+1-j),y(i,end+1-j),u(i,end+1-j),...
                                                                  v(i,end+1-j),x(i+1,j-1), y(i+1,j-1),u(i+1,j-1),...
                                                                         v(i+1,j-1),To,R,gamma,del,Po) ;% Generates Int points.
             
              
             modv(i+1,j) = (u(i+1,j).^2 + v(i+1,j).^2).^0.5 ;            % Speed at nodal point.
             al(i+1,j) = (ao.^2 - ((gamma - 1)/2)*modv(i+1,j).^2).^0.5 ; % Relative acoustic speed at nodal location.
             M(i+1,j) = modv(i+1,j)./al(i+1,j) ; 
             
             Pm(i+1,j) = ((1 + f*M(i+1,j).^2).^(-gamma/(2*f)))*Po ;
             Tm(i+1,j) = ((1 + f*M(i+1,j).^2).^(-1))*To ;
             Rho(i+1,j) = Pm(i+1,j)/(R*Tm(i+1,j)) ;
             
             k(i+1,j) = (Rho(i+1,j)*modv(i+1,j)*sind(atand(((y(i+1,j)-y(i+1,j-1))/(x(i+1,j)-x(i+1,j-1)))) - ...
                                                      (atand(v(i+1,j)/u(i+1,j))))*y(i+1,j))/(sind(atand(((y(i+1,j)...
                                                                                  -y(i+1,j-1))/(x(i+1,j)-x(i+1,j-1)))))) ;
              
             mdot(i+1-p,j-1) = pi*(k(i+1,j)+k(i+1,j-1))*(y(i+1,j)-y(i+1,j-1)) ;
             
             Mdot = sum(mdot(i+1-p,:)) + sum(mline(1:i-p+1)) ; 
             if Mdot > Mdotnumeric
                mc = sum(mdot(i+1-p,1:end-1)) + sum(mline(1:i-p+1)) ;
                xt = x(i+1,j-1) ;
                yt = y(i+1,j-1) ;
                ut = u(i+1,j-1) ; 
                vt = v(i+1,j-1) ;
            for w = 1:25
                 xmid = (xt+x(i+1,j))/2 ;
                 ymid= (yt+y(i+1,j))/2 ;
                 umid = (ut+u(i+1,j))/2 ;
                 vmid = (vt+v(i+1,j))/2 ;
                 
                 modva = (umid^2 + vmid^2)^0.5 ;
                 ala = (ao.^2 - ((gamma - 1)/2)*modva.^2).^0.5 ;
                 ma = modva/ala ;
                 Pma = ((1 + f*ma.^2).^(-gamma/(2*f)))*Po ;
                 Tma = ((1 + f*ma.^2).^(-1))*To ;
                 Rhoa = Pma/(R*Tma) ;
     
                 ka = (Rhoa*modva*sind(atand(((ymid-yt)/(xmid-xt))) - ...
                              (atand(vmid/umid)))*ymid)/(sind(atand(((ymid-yt)/(xmid-xt))))) ;
     
                 mdota = pi*(k(i+1,j-1)+ka)*(ymid-yt) ;
                
                  if (mc + mdota) > Mdotnumeric
                        x(i+1,j) = xmid ;
                        y(i+1,j) = ymid ;
                        u(i+1,j) = umid ;
                        v(i+1,j) = vmid ;
                        
                        modv(i+1,j) = modva ;
                        M(i+1,j) = ma ;
                        Pm(i+1,j) = Pma ;
                        Tm(i+1,j) =  Tma ;
                        Rho(i+1,j) = Rhoa ;
                       
                        mdot(i+1-p,end) = mdota ;
                  elseif(mc + mdota) < Mdotnumeric
                        xt = xmid ;
                        yt = ymid ;
                        ut = umid ;
                        vt = vmid ;
                        mc = mc + mdota ;
                  end
            end
             break % For breaking the j loop using the mass flow rate condition.
             end
         end
     end
 end
     for i = p+1:p+n-1
       for j = 2:q
          
            [x(i+1,j),y(i+1,j),u(i+1,j),v(i+1,j)] =   MOC_int_point(x(i,j),y(i,j),u(i,j),v(i,j),...
                                                     x(i+1,j-1),y(i+1,j-1),u(i+1,j-1),v(i+1,j-1),To,R,gamma,del,Po) ;
                                                 
                                                 
             modv(i+1,j) = (u(i+1,j).^2 + v(i+1,j).^2).^0.5 ;            % Speed at nodal point.
             al(i+1,j) = (ao.^2 - ((gamma - 1)/2)*modv(i+1,j).^2).^0.5 ; % Relative acoustic speed at nodal location.
             M(i+1,j) = modv(i+1,j)./al(i+1,j) ; 
             
             Pm(i+1,j) = ((1 + f*M(i+1,j).^2).^(-gamma/(2*f)))*Po ;
             Tm(i+1,j) = ((1 + f*M(i+1,j).^2).^(-1))*To ;
             Rho(i+1,j) = Pm(i+1,j)/(R*Tm(i+1,j)) ;
             
             k(i+1,j) = (Rho(i+1,j)*modv(i+1,j)*sind(atand(((y(i+1,j)-y(i+1,j-1))/(x(i+1,j)-x(i+1,j-1)))) - ...
                                                (atand(v(i+1,j)/u(i+1,j))))*y(i+1,j))/(sind(atand(((y(i+1,j)-y(i+1,j-1))...
                                                                                                   /(x(i+1,j)-x(i+1,j-1)))))) ;
              
          
              mdot(i+1-p,j-1) = pi*(k(i+1,j)+k(i+1,j-1))*(y(i+1,j)-y(i+1,j-1)) ;
             
             Mdot = sum(mdot(i+1-p,:)) + sum(mline(1:i-p+1)) ; 
             
             if Mdot > Mdotnumeric
                 z = nnz(mdot(i+1-p,:));
                mc = sum(mdot(i+1-p,1:z-1)) + sum(mline(1:i-p+1))  ;
            % While abs(Mdot - Mdotnumeric) > 10^-1
            xt = x(i+1,j-1) ;
            yt = y(i+1,j-1) ;
            ut = u(i+1,j-1) ; 
            vt = v(i+1,j-1) ;
            
            for w = 1:25
                 xmid = (xt+x(i+1,j))/2 ;
                 ymid= (yt+y(i+1,j))/2 ;
                 umid = (ut+u(i+1,j))/2 ;
                 vmid = (vt+v(i+1,j))/2 ;
                 
                 modva = (umid^2 + vmid^2)^0.5 ;
                 ala = (ao.^2 - ((gamma - 1)/2)*modva.^2).^0.5 ;
                 ma = modva/ala ;
                 Pma = ((1 + f*ma.^2).^(-gamma/(2*f)))*Po ;
                 Tma = ((1 + f*ma.^2).^(-1))*To ;
                 Rhoa = Pma/(R*Tma) ;
     
                 ka = (Rhoa*modva*sind(atand(((ymid-yt)/(xmid-xt))) - ...
                              (atand(vmid/umid)))*ymid)/(sind(atand(((ymid-yt)/(xmid-xt))))) ;
     
                 mdota = pi*(k(i+1,j-1)+ka)*(ymid-yt) ;
                
                  if (mc + mdota) > Mdotnumeric
                        x(i+1,j) = xmid ;
                        y(i+1,j) = ymid ;
                        u(i+1,j) = umid ;
                        v(i+1,j) = vmid ;
                     
                        modv(i+1,j) = modva ;
                        M(i+1,j) = ma ;
                        Pm(i+1,j) = Pma ;
                        Tm(i+1,j) =  Tma ;
                        Rho(i+1,j) = Rhoa ;
                        
                        mdot(i+1-p,end+1-z) = mdota ;
                  elseif(mc + mdota) < Mdotnumeric
                        xt = xmid ;
                        yt = ymid ;
                        ut = umid ;
                        vt = vmid ;
                        mc = mc + mdota ;
                  end
            end
             break % For breaking the j loop using the mass flow rate condition.
             end
       end
     end
     

%% Plots the mesh of the turning section.

 hold on  
 for i = p+1:p+n
     x1(i-p) = x(i,nnz(x(i,:)));
     y1(i-p) = y(i,nnz(y(i,:)));
 end
 plot(x1,y1,'-k','linewidth',2);
 hold on % Here we plot the internal mesh of the Turning section.
 for i = p+1:p+n-1 
      for j = 2:nnz(x(i+1,:))
              hold on 
               plot([x(i,j),x(i+1,j)],[y(i,j),y(i+1,j)],'k-','linewidth',2)
              hold on
               plot([x(i+1,j-1),x(i+1,j)],[y(i+1,j-1),y(i+1,j)],'k-','linewidth',2)
      end
 end
%}
%% This joins the circular arc extension to the turning contour.
 
 hold on 
 for i = 1:nnz(x(p+1,:))
      xj(i) = x(p,nnz(x(p,:))-i+1);
      yj(i) = y(p,nnz(y(p,:))-i+1);
      xk(i) = x(p+1,i);
      yk(i) = y(p+1,i);
 end
 for i = 1:nnz(xk)
     plot([xj(i),xk(i)],[yj(i),yk(i)],'linewidth',2)
 end
for i = 1:nnz(xk)-1
    plot([xk(i+1),xk(i)],[yk(i+1),yk(i)],'k-','linewidth',2)
end
hold on 
plot([x(p),x(p+1,nnz(x(p+1,:)))],[y(p),y(p+1,nnz(x(p+1,:)))],'k-','linewidth',2)
%} 
%% Here we plot the nozzle contour only and no mesh.  

figure(3)
clear x1 
clear y1
  for i = n1+1:a
    x1(i-n1) = x(i,1);
    y1(i-n1) = y(i,1);
  end
  hold on 
  plot(x1,y1,'k-');
  clear x1
  clear y1
   for i = p+1:p+n
     x1(i-p) = x(i,nnz(x(i,:)));
     y1(i-p) = y(i,nnz(y(i,:)));
   end
  plot(x1,y1,'*');
 hold on 
 plot(x1,y1,'k-');
 hold on 
 plot([x(p),x(p+1,nnz(x(p+1,:)))],[y(p),y(p+1,nnz(x(p+1,:)))],'k-','linewidth',2)
 xlim([0 xivlbar(1)+ 15]);
 ylim([0 th+2.0]);
 title(sprintf('Nozzle contour for Mach=%d and cp/cv =%d ',mi,gamma));
 xlabel('Axial length');
 ylabel('Radial length');
 
%} 
 
%% Variation of Temperature of Wall VS the Nozzle length.
%{
 figure(4)
 clear x1 
 clear t1
   for i = n1+1:a
     x1(i-n1) = x(i,1);
     t1(i-n1) = Tm(i,1);
   end
   for i = p+1:p+n
      x1(i-n1) = x(i,nnz(x(i,:)));
      t1(i-n1) = Tm(i,nnz(y(i,:)));
   end
  plot(x1,t1,'*');
  title(sprintf('Variation of Temperature along the Nozzle-wall for Mach=%d ',mi));
  xlabel('Axial length of the Nozzle');
  ylabel('Temperature(K)');
%}

%% Variation of Mach Number along the Nozzle wall against the Nozzle Axis.
%{
 figure(5)
 clear x1 
 clear m1
   for i = n1+1:a
     x1(i-n1) = x(i,1);
     m1(i-n1) = M(i,1);
   end
   for i = p+1:p+n
      x1(i-n1) = x(i,nnz(x(i,:)));
      m1(i-n1) = M(i,nnz(y(i,:)));
   end
  plot(x1,m1,'.');
  title(sprintf('Variation of Mach number along the Nozzle-wall for Mach=%d ',mi));
  xlabel('Axial length of the Nozzle');
  ylabel('Mach-Number');
  xlim([-0.1 15]) ;
  hold on 
 
  % The Axis-Symmetry points. Centreline mach number
  clear x1
  clear m1
 %  figure(6)
  for i = 1:n1
        x1(i)= x(i,1);
        m1(i) = M(i,1);
  end
  plot(x1,m1,'o')
    hold on
    clear x2
    clear m2
  for i = n1+1:a
      x2(i-n1) = x(i,nnz(x(i,:)));
      m2(i-n1) = M(i,nnz(M(i,:)));
   end
  plot(x2,m2,'*');
 hold on 
 plot([x1(end) x2(1)],[m1(end) m2(1)]);
  clear x1 
  clear m1
  for i = a+1:a+n
      x1(i-a) = x(i,1);
      m1(i-a) = M(i,1);
   end
  plot(x1,m1,'*');
%}  
%% Variation of Pressure along the Nozzle wall against Nozzle length.
%{

 clear x1 
 clear p1
   for i = n1+1:a
     x1(i-n1) = x(i,1);
     p1(i-n1) = Pm(i,1)/Po;
   end
   for i = p+1:p+n
      x1(i-n1) = x(i,nnz(x(i,:)));
      p1(i-n1) = Pm(i,nnz(y(i,:)))/Po;
   end
  plot(x1,p1,'.');

%}

%% Variation of pressure along symmetry axis vs  nozzle length.
%{
figure(7)
  hold on 
  clear x1 
  clear p  
 for i = 1:n1 
 x1(i)=x(i,1);
 p(i) = Pm(i,1)/Po ;
 end
 plot(x1,p,'.-')
 hold on 
 clear x2
 clear p2 
 for i = n1+1:a
    x2(i-n1) = x(i,nnz(x(i,:)));
    p2(i-n1) = Pm(i,nnz(x(i,:)))/Po;
 end
 plot(x2,p2,'.-')
 hold on 
 plot([x1(end) x2(1)],[p(end) p2(1)]);
  title(sprintf('Variation of Pressure ratio along the Nozzle-wall and axis of symmetry for Mach=%d ',mi));
  xlabel('Axial length of the Nozzle');
  ylabel('Pressure Ratio');
  xlim([0 8]);

%}
 %% This code is used for creating Iso-Mach line, manually input two ranges of Mach
              %that are close to each other.
%{
  
  hold on
  clear r 
  clear c 
  clear xp 
  clear yp  
  [r,c] = find(M>1.496 & M<1.504);
    for i = 1:length(r) 
        xp(i) = [x(r(i),c(i))];
        yp(i) = [y(r(i),c(i))];
    end
 plot(xp,yp,'.')

              
%}   
%% Here we plot the two y axis-one mach number and other pressure ratio vs length of nozzle.
%{
 figure
 clear x1 
 clear m1
 
   for i = n1+1:a
     x1(i-n1) = x(i,1);
     m1(i-n1) = M(i,1);
   end
   for i = p+1:p+n
      x1(i-n1) = x(i,nnz(x(i,:)));
      m1(i-n1) = M(i,nnz(y(i,:)));
   end
   
  clear x1 
  clear p1
   for i = n1+1:a
     x1(i-n1) = x(i,1);
     p1(i-n1) = Pm(i,1)/Po;
   end
   for i = p+1:p+n
      x1(i-n1) = x(i,nnz(x(i,:)));
      p1(i-n1) = Pm(i,nnz(y(i,:)))/Po;
   end
   [AX,H1,H2] = plotyy(x1,m1,x1,p1) ;
   set(get(AX(1),'Ylabel'),'string','Mach number','Fontsize',14)
   set(get(AX(2),'Ylabel'),'string','P/P_o','Fontsize',14)
%}

  


         
         
      
 
 
                 
                 
             





