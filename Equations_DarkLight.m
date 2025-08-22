%--------------------------------------------------------------------------
% LICENSE
% Copyright (C) 2025  Davide Moia
% Max Planck Institute for Solid State Research
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%--------------------------------------------------------------------------
function F = Equations_DarkLight(x, par, Rec_Model, Redox_Model)

pI2 = par(1);
kf_n_i = par(2);
kf_p_i = par(3);
kb_n_i = par(4);
kb_p_i = par(5);
kf_n_v = par(6);
kf_p_v = par(7);
kb_n_v = par(8);
kb_p_v = par(9);
kf_sg_i = par(10);
kf_sg_v = par(11);
kb_sg_i = par(12);
kb_sg_v = par(13);
KaF = par(14);
kion_rec = par(15);
ni = par(16);
n1 = par(17);
p1 = par(18);
taun = par(19);
taup = par(20);
Gext = par(21);
krad = par(22);
Cn = par(23);
Cp = par(24);

Ki = kf_sg_i/kb_sg_i;
Kv = kf_sg_v/kb_sg_v;

rad_coeff = Rec_Model(1);
nr_coeff = Rec_Model(2);
aug_coeff = Rec_Model(3);

i_coeff = Redox_Model(1);
v_coeff = Redox_Model(2);

%Electroneutrality
F(1) = x(1) + x(3) - x(4) - x(2);

if (v_coeff~=i_coeff) %Assumes i_coeff = 1
    %Iodide vacancy balance
    F(2) = x(3)*x(4)-KaF;
    %Iodide interstitial balance
    F(3) = -kf_p_i*x(1)*x(4) + kb_p_i*(pI2^0.5)/Ki - kf_n_i*x(4) + kb_n_i*x(2)*(pI2^0.5)/Ki;
    %Iix
    F(4) = x(6) - (pI2^0.5)/Ki;
    %VIx
    F(5) = x(5) - Kv/(pI2^0.5);
else
    %Iodide vacancy balance
    F(2) = kf_n_v*x(5) - kb_n_v*x(3)*x(2) + kf_p_v*x(5)*x(1) - kb_p_v*x(3) - kion_rec*(x(3)*x(4)-KaF);
    %Iodide interstitial balance
    F(3) = -kf_p_i*x(1)*x(4) + kb_p_i*x(6) - kf_n_i*x(4) + kb_n_i*x(2)*x(6) - kion_rec*(x(3)*x(4)-KaF);
    %Iix balance
    F(4) = kf_p_i*x(1)*x(4) - kb_p_i*x(6) + kf_n_i*x(4) - kb_n_i*x(2)*x(6) - kf_sg_i*x(6) + kb_sg_i*(pI2^0.5);
    %VIx balance
    F(5) = -(kf_p_v*x(1) + kf_n_v)*x(5) + (kb_p_v + kb_n_v*x(2))*x(3) + kf_sg_v - kb_sg_v*x(5)*(pI2^0.5);
end
% Electronic disorder
% Rrad, RSRH, RAug, R from redox with Iodine defects
UI_i = i_coeff*(kb_n_i*x(2)*x(6) - kf_n_i*x(4));
UI_v = v_coeff*(kb_n_v*x(3)*x(2) - kf_n_v*x(5));
F(6) = Gext - (x(2)*x(1)-ni^2)*(rad_coeff*krad + nr_coeff/((taup)*(x(2)+n1)+(taun)*(x(1)+p1)) + aug_coeff*(Cp*x(1) + Cn*x(2))) - UI_i - UI_v;

end
