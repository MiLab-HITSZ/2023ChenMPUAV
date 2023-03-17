function np=getNoiseimpact(h)
d=88;
Lh=55;
np = 1.5499e+04*Lh./((h./0.348).^2+d^2);
if np<40
    np=0;
end
end