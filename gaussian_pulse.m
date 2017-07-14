function pulse = gaussian_pulse(t,t_mt,amp,bw)
% gaussian_pulse calculates gaussian pulse of t_mt s duration

% inputs: 
%   t - time, s
%   t_mt - duration of pulse, s
%   amp - amplitude of pulse, s-1
%   bw - bandwidth of pulse, Hz
%   options - additional parameters, including options.bw - bandwidth of 
%             the pulse
% output: 
%   pulse - computed pulse shape

% sigma_time_domain*sigma_freq_domain = 1/(2*pi)
% bandwidth defined as the half power point: for an arbitrary cut-off 
% value 1/c for magnitude of pulse the cut-off frequency is given 
% by BW/2 = sqrt(2*ln(c))*sigma_freq_domain
% therefore sigma_time_domain = sqrt(2*ln(2))/(pi*BW)

sigma = sqrt(2*log(2)) ./ (pi*bw);
pulse = amp' .* exp( -((t-(t_mt/2)).^2)./(2*sigma.^2));
pulse((t < 0 | t>t_mt)) = 0;

return

