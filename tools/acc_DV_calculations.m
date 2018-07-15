Acor(:,i+1) = 2*fast_cross(omega_TAR_INT(:,i+1),rdot_SVC_TAR(:,i+1));
A_COR(i+1) = norm(Acor(:,i+1));
Aang(:,i+1) = fast_cross(omega_dot_TAR_INT(:,i+1),r_SVC_TAR(:,i+1));
A_ANG(i+1) = norm(Aang(:,i+1));
Acen(:,i+1) = fast_cross(omega_TAR_INT(:,i+1), fast_cross(omega_TAR_INT(:,i+1),r_SVC_TAR(:,i+1)));
A_CEN(i+1) = norm(Acen(:,i+1));
% if i == 1
%     Alin(:,i+1) = (rdot_SVC_TAR(:,i+1) - rdot_SVC_TAR(:,i)) /dt;
% elseif i == 2
%     Alin(:,i+1) = (rdot_SVC_TAR(:,i+1) - rdot_SVC_TAR(:,i-1)) /(2*dt);
% else
%     Alin(:,i+1) = (rdot_SVC_TAR(:,i+1) - rdot_SVC_TAR(:,i-2)) /(3*dt);
% end
Alin(:,i+1) = (rdot_SVC_TAR(:,i+1) - rdot_SVC_TAR(:,i)) /dt;
A_LIN(i+1) = norm(Alin(:,i+1));

Atot(:,i+1) = Alin(:,i+1) + Acor(:,i+1) + Aang(:,i+1) + Acen(:,i+1);
A_TOT(i+1) = sqrt(Atot(1,i+1)^2 + Atot(2,i+1)^2 + Atot(3,i+1)^2);

A_LIN_RATIO(i+1) = A_LIN(i+1)/A_TOT(i+1);
A_COR_RATIO(i+1) = A_COR(i+1)/A_TOT(i+1);
A_ANG_RATIO(i+1) = A_ANG(i+1)/A_TOT(i+1);
A_CEN_RATIO(i+1) = A_CEN(i+1)/A_TOT(i+1);

DV_angular(i+1)     = DV_angular(i) + A_ANG(:,i+1)*dt;
DV_centripetal(i+1) = DV_centripetal(i) + A_CEN(i+1)*dt;
DV_linear(i+1)      = DV_linear(i) + A_LIN(:,i+1)*dt;
DV_coriolis(i+1)    = DV_coriolis(i) + A_COR(i+1)*dt;
DV_total(i+1) = DV_total(i) + A_TOT(i+1)*dt;