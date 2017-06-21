function result = lorentzian(delta,T2)
% lorentzian lineshape function
% inputs: 
%   T2 - relaxation time of restiricted pool in s (scalar value)
%   delta - linear frequency offset of RF pulse, Hz (scalar value/vector)
% output: 
%   calculate lineshape (scalar value/vector)

result = T2./(pi*(1+(2*pi*delta*T2).^2));

