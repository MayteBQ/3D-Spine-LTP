function S_int = solve_system_RK_3D_iso(S,ff,actin_dyn,delta_t,P)
% solve system using RK4
    if actin_dyn == 0 
        k1 = delta_t*P.zeta*(membrane_force_3D_iso(S,ff,P));
        k2 = delta_t*P.zeta*(membrane_force_3D_iso(S+k1./2,ff,P));
        k3 = delta_t*P.zeta*(membrane_force_3D_iso(S+k2./2,ff,P));
        k4 = delta_t*P.zeta*(membrane_force_3D_iso(S+k3,ff,P));
    else
        k1 = delta_t*P.zeta*(membrane_force_3D_iso(S,ff,P)+f_fil_3D(S,P));
        k2 = delta_t*P.zeta*(membrane_force_3D_iso(S+k1./2,ff,P)+f_fil_3D(S+k1./2,P));
        k3 = delta_t*P.zeta*(membrane_force_3D_iso(S+k2./2,ff,P)+f_fil_3D(S+k2./2,P));
        k4 = delta_t*P.zeta*(membrane_force_3D_iso(S+k3,ff,P)+f_fil_3D(S+k3,P));
    end
      
    S_int = S + (k1+2*k2+2*k3+k4)/6;

end