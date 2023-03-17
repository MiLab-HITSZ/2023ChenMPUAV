function cos=getcos2(h,delta_h,W,t)
delta_h=max(delta_h,0);
area_blade = 0.1;
n=4;
g=9.8;
rou=htorou(h);
cos=max(W^(3/2)*sqrt(g^3/2/rou/area_blade/n)*t+delta_h*W,0);
end
