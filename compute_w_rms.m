function w_rms = compute_w_rms(pulse)
% Approximate MT pulse using Yarnykh and Yuan?s model, replacing MT shaped
% pulse with rectangular pulse of the same duration and equivalent 
% average power

% input: pulse - set of pulse parameters
% output: averaged omega1 value

t_mt = pulse.t_mt;
omega1 = pulse.omega1;
integral_value = integral(omega1.^2, 0, t_mt);
w_rms = sqrt(integral_value / t_mt);

end

