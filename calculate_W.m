function W = calculate_W(delta,T2r,lineshape, pulse)
% calculate mean saturation rate
% input:
%   delta - linear frequency offset of RF pulse, Hz (scalar value/vector)
%   T2 - relaxation time of restiricted pool in s (scalar value)
%   lineshape - assumed lineshape of restricted pool
%   pulse - data structure describing MT pulse
% output: 
%   W - mean saturation rate for given parameters
G = 0;
switch lineshape
    case 'gaussian'
        G = gaussian(delta, T2r);
    case 'lorentzian'
        G = lorentzian(delta, T2r);
    case 'superlorentzian'
        G = superlorentzian(delta, T2r);
    case 'superlorentzian_res'
        G = superlorentzian_res(delta, T2r);
end


end

W = zeros(length(Pulse),1);

for ii = 1:length(Pulse)
    omega2 = Pulse(ii).omega2;
    Trf = Pulse(ii).Trf;
    int = integral(omega2, 0, Trf, 'ArrayValued', true);
    W(ii) = G * pi/Trf * int;
end