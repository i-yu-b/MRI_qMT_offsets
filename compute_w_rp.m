function w_rp = compute_w_rp(TR, pulse)
% Approximate MT pulse using Sled and Pike RP approach
% w_rp function computes cw pulse of pulse duration equals to 
% full-width-at-half-maximum of the shaped MT pulse in pulse sequence of
% duration equals to t_mt

% input: TR - repetition time, pulse - set of pulse parameters
% output: averaged omega1 value

t_mt = pulse.t_mt;
omega1 = pulse.omega1;
integral_value = integral(omega1.^2, 0, t_mt);

if strcmp(pulse.shape,'hard')
    tau = t_mt;
else
    t = 0:2000:t_mt;
    tau = fwhm(t,omega1);
end

w_rp = sqrt(integral_value/tau);

end
