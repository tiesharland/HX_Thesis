function [T_mean_h, T_mean_c, T_c_i, T_c_o, T_h_i, T_h_o, C_h, C_c, C_star] = HX_deltaT (T3, M_dot_3,M_dot_coolant_max,T_out_fc, T_in_fc)

%------Hot side - Coolant loop; Cold side - Air flow ; This will be the convention used here after-------%

Cp_air_hx_inlet = 1005 + 0.1*(T3-288) - 0.0002*(T3-288).^2; %calculation of cp of air at ambient temperature
effectiveness = 0.8; %assumption based on HEX design procedures
C_c = M_dot_3*Cp_air_hx_inlet;
C_h = M_dot_coolant_max*cp_EG_50_50(T_out_fc);
C_star = C_c/C_h;

T_c_i = T3; %K, cool side i.e air side inlet temp into HX is outlet temperature of diffuser
T_h_i = T_out_fc;%K, hot side i.e. coolant side inlet temp into HX is fuel cell outlet temp.

T_h_o = T_in_fc; %K, hot side i.e coolant side outlet temp from HX
T_c_o = T_c_i+(effectiveness*(C_c*(T_h_i-T_c_i)/C_c)); %K, cool side i.e air side outlet temp from HX

T_mean_h = (T_h_i+T_h_o)/2;
T_lm = ((T_mean_h-T_c_o)-(T_mean_h-T_c_i))/log((T_mean_h-T_c_o)/(T_mean_h-T_c_i));
T_mean_c = T_mean_h-T_lm;


end
