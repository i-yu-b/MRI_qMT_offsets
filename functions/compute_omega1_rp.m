function omega1_rp = compute_omega1_rp(pulse)
%{
omega1_rp = compute_omega1_rp(pulse) function computes mean omega1 value 
of shaped  MT pulse using Sled and Pike RP approach: MT shaped pulse of 
duration t_mt  approximated by cw rectangular pulse of duration equals to 
full-width-at-half-maximum of the shaped MT pulse and equivalent average 
power

input:
    pulse - data structure containig all MT pulse properties
output:
    omega1_rp - averaged omega1 value, [rad/s]
%}

TR = pulse.TR;
t_mt = pulse.t_mt;
omega1_sq_func = pulse.omega1_sq_func;  % function handle describing power 
                                        % of MT pulse depending on time
omega1_func = pulse.omega1_func;  % function handle describing amplitude of 
                                  % MT pulse depending on time
                                        
integral_value = integral(omega1_sq_func, 0, t_mt);

if strcmp(pulse.shape,'<hard>')
    tau = t_mt;
else
    t = linspace(0, t_mt, 1000);
    tau = fwhm(t,omega1_func(t)) % find full-width-at-half-maximum value 
end

omega1_rp = sqrt(integral_value/tau) * pulse.amp * 2 *pi;

end
