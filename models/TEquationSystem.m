function F = TEquationSystem(K1, K2, K3, K4, Ku1, Ku2, Q10, Q20, T_inf, x)
    T10 = x(1);
    T20 = x(2);
    
    F(1) = K1*T_inf+K2*T_inf^4+Ku1*Q10 + T10*(-K1-K3) + T20*K3 + T10^4*(-K2-K4) + T20^4*K4;
    F(2) = K1*T_inf+K2*T_inf^4+Ku2*Q20 + T20*(-K1-K3) + T10*K3 + T20^4*(-K2-K4) + T10^4*K4;
end