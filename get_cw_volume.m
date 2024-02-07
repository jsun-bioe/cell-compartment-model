function Vb = get_cw_volume(Vtot, Vbulge, kLpp, kCW, Vr, V1i)
%GET_CW_VOLUME 
%   

Vb = (kLpp*(Vtot-Vbulge)/V1i - kLpp + kCW) / (kCW/Vr + kLpp/V1i);

end

