%--------------------------------------------------------------------------
% LICENSE
% Copyright (C) 2025  Davide Moia
% Max Planck Institute for Solid State Research
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%--------------------------------------------------------------------------
function F = Equations_DarkLight_sgeq(x, par, Rec_Model, Redox_Model)

pI2 = par(1);
kf_n_i = par(2);
kf_p_i = par(3);
kb_n_i = par(4);
kb_p_i = par(5);
kf_n_v = par(6);
kf_p_v = par(7);
kb_n_v = par(8);
kb_p_v = par(9);
KaF = par(10);
kion_rec = par(11);
ni = par(12);
n1 = par(13);
p1 = par(14);
taun = par(15);
taup = par(16);
Gext = par(17);
krad = par(18);
Cn = par(19);
Cp = par(20);
Ki = par(21);
Kv = par(22);

rad_coeff = Rec_Model(1);
nr_coeff = Rec_Model(2);
aug_coeff = Rec_Model(3);

i_coeff = Redox_Model(1);
v_coeff = Redox_Model(2);

%Electroneutrality
F(1) = x(1) + x(3) - x(4) - x(2);

%If only redox involving interstitials are active the rate equation for
%interstitials is used, neglecting the thermal generation and recombination
%of anti-Frnkel pairs (they are at equilibrium)
if i_coeff ~= v_coeff
    F(2) = -kf_p_i*x(1)*x(4) + kb_p_i*(pI2^0.5)/Ki - kf_n_i*x(4) + kb_n_i*x(2)*(pI2^0.5)/Ki;
%If redox involving both interstitials and vacancies are active, the
%complete rate equation for interstitials is used for F(2), including 
%thermal generation and recombination
else
    F(2) = -kf_p_i*x(1)*x(4) + kb_p_i*(pI2^0.5)/Ki - kf_n_i*x(4) + kb_n_i*x(2)*(pI2^0.5)/Ki - kion_rec*(x(3)*x(4)-KaF);
end

%If only redox involving interstitials are active, F(3) sets the anti-Frenkel equilibrium 
if i_coeff ~= v_coeff
    F(3) = x(3)*x(4)-KaF;
%If redox involving both interstitials and vacancies are active, the
%complete rate equation for vacancies is used fro F(3)
else
    F(3) = v_coeff*((kf_p_v*x(1)*Kv + kf_n_v*Kv)/pI2^0.5 - (kb_p_v + kb_n_v*x(2))*x(3) - kion_rec*(x(3)*x(4)-KaF));
end

% Electronic disorder
% Rrad, RSRH, RAug, R from redox with Iodine defects
F(4) = Gext - (x(2)*x(1)-ni^2)*(rad_coeff*krad + nr_coeff/((taup)*(x(2)+n1)+(taun)*(x(1)+p1)) + aug_coeff*(Cp*x(1) + Cn*x(2))) - (i_coeff*(kb_n_i*x(2)*pI2^0.5/Ki - kf_n_i*x(4)) + v_coeff*(kb_n_v*x(3)*x(2) - kf_n_v*Kv*pI2^-0.5));

end
