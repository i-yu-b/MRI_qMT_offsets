function pulse_shape = gaussian_pulse(t,t_mt, bw)
%{
pulse_shape = hard_pulse(t,t_mt,amp) function provides function handle to
describe lineshape(t) of hard pulse of t_mt s duration and amplitude = 1

% inputs: 
%   t - time, s
%   t_mt - duration of pulse, s
%   bw - bandwidth of pulse, Hz
% output: 
%   pulse_shape - function handle of computed pulse shape 
%}

% sigma_time_domain*sigma_freq_domain = 1/(2*pi)
% bandwidth defined as the half power point: for an arbitrary cut-off 
% value 1/c for magnitude of pulse the cut-off frequency is given 
% by BW/2 = sqrt(2*ln(c))*sigma_freq_domain
% therefore sigma_time_domain = sqrt(2*ln(2))/(pi*BW)

sigma = sqrt(2*log(2)) ./ (pi*bw);
pulse_shape = exp( -((t-(t_mt/2)).^2)./(2*sigma.^2));
pulse_shape((t < 0 | t>t_mt)) = 0;

return

