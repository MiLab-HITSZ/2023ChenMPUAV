function V = getV(h,data)
p1=2*data.m*data.g/data.R_I/data.S_hit/data.rou_a;
p2=(1-exp(-(h*data.R_I*data.S_hit*data.rou_a)/data.m));
V = sqrt(p1*p2);
end