function result = gaussian(delta,T2)
% gaussian lineshape function
% inputs: 
%   T2 - relaxation time of restiricted pool in s (scalar value)
%   delta - linear frequency offset of RF pulse, Hz (scalar value/vector)
% output: 
%   calculate lineshape (scalar value/vector)

result = T2/sqrt(2*pi)*exp(-(2*pi*delta*T2).^2/2);
