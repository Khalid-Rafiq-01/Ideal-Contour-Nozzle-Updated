function[xc,yc] = circle_arc(rc,theta,th) 
xc = rc*sind(theta) ;
yc = th + rc*(1-cosd(theta)) ;
end