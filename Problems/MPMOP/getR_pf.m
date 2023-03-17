function R_pf = getR_pf(V,data)
E_imp=data.m*V*V/2;
R_pf=1/(1+sqrt(data.alpha/data.beta)*(data.beta/E_imp)^(1/(4*data.S_c)));
end