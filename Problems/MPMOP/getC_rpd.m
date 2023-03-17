function C_rpd = getC_rpd(h,data)
h = max(h,exp(data.miu));
C_rpd=1./(h.*data.sigma.*sqrt(2.*pi)).*exp(-(log(h)-data.miu).^2./(2*data.sigma^2));
end