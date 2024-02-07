function Vbulge = get_bulge_volume(Vtot, V1i, Vcw, Emax, blpp)
%GET_BULGE_VOLUME 
%   

elpp = (Vtot-Vcw)/V1i - 1;

s = exp(-blpp*(elpp-Emax))/(1+exp(-blpp*(elpp-Emax)));

Vperi = V1i * (1 + elpp*s);

Vbulge = Vtot - Vcw - Vperi;

end

