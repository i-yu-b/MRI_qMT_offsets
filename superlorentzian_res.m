function result = superlorentzian_res(delta,T2)
% superlorentzian lineshape function in case of the presence on-resonance
% e.g. delta = 0
% cutoff frequency for on-resonance case is 1.5 kHz
% inputs: 
%   T2 - relaxation time of restiricted pool in s (scalar value)
%   delta - linear frequency offset of RF pulse, Hz (scalar value/vector)
% output: 
%   calculate lineshape (scalar value/vector)

    % set up cutoff frequency at 1.5 kHz
    cutoff = 1.5e3;
    % finding cutoff for negative values
    idx_plus = find(delta >= cutoff, 1, 'first'); 
    % finding cutoff for positive values
    idx_minus = find(delta <= -cutoff, 1, 'last'); 
    
    % function for calculating superlorentzian linewidth
    integrand = @(u,delta,T2) sqrt(2/pi)*T2./...
              abs(3*u.^2-1).*exp(-2*((2*pi*delta*T2)'*(1./(3*u.^2-1))).^2);

    if (~isempty(idx_plus) && ~isempty(idx_minus)) % both cutoffs found 
      
        % exclude cutoff interval from delta values for calculations
         if isrow(delta) %delta is a row vector
             delta_calc = [delta(1:idx_minus) delta(idx_plus:end)];
         else %delta is a column vector
             delta_calc = [delta(1:idx_minus); delta(idx_plus:end)]; 
         end
         
    elseif (~isempty(idx_plus)) % cutoff for positive values found 
        
        % exclude cutoff interval from delta values for calculations
        delta_calc = delta(idx_plus:end);
        
    elseif (~isempty(idx_minus)) % cutoff for negative values found
        
        % exclude cutoff interval from delta values for calculations
        delta_calc = delta(1:idx_minus);
      
    else % no on-resonanve values presented
        
        % execute superlorentzian function and end current function
        result = superlorentzian(delta,T2);
        return
        
    end
          
    % calculate integral in off-resonance range
    result_offres = integral(@(u)integrand(u,delta_calc,T2),0,1,...
                             'ArrayValued',true); 
    % interpolate result to cutoff interval
    result = interp1(delta_calc,result_offres,delta,'spline');
  
    
    