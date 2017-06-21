function R1f = compute_R1(parameters, R1_obs)

kfr  = parameters.kfr;
F  = parameters.F;
R1r = parameters.R1r;

R1f = R1_obs - kfr*(R1r - R1_obs) / (R1r - R1_obs + kfr/F);

end