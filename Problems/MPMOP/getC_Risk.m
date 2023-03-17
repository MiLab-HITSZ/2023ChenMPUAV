function C_Risk = getC_Risk(Risk,Density,data)
C_Risk = data.P_crash*Density*Risk*data.S_hit;
end