function fwhm_value = fwhm(x,y)
% fwhm function calculates full-width at half-maximum (FWHM) of the 
% signal described by y(x)
% input: x, y - 
% output: fwhm_value - FWHM of y(x) in units of x

% normalize the signal
y = y/max(y);

% use min-max scaler if necessary (all data point are above 0.5 level)
if isempty(find(y<0.5)) 
    y = (y-min(y))/ (max(y)-min(y));
end

n = length(y);
plot(x, y, '-ok');
hold on;

% find indexs of left and right points corresponding to half-maximum.
indexes = find(diff(sign(y-0.5))~=0,2,'first');
idx_left = indexes(1)+1;
idx_right = indexes(2)+1;


if ~isempty(idx_left)
   % interpolate the results
    coef_left = (y(idx_left)-0.5)/(y(idx_left)-y(idx_left-1));
    x_left = x(idx_left-1)+coef_left*(x(idx_left)-x(idx_left-1));
else 
    fwhm_value=NaN;
    return
end

 
if ~isempty(idx_right)
   % interpolate the results
    coef_right = (0.5 - y(idx_right))/(y(idx_right-1)-y(idx_right));
    x_right = x(idx_right-1)+coef_right*(x(idx_right)-x(idx_right-1));
    y_right = y(idx_right-1)-coef_right*(x_right-x(idx_right-1));
else 
    fwhm_value = NaN;
    return
end

fwhm_value = x_right - x_left;

plot([x_left, x_right], [0.5,0.5],'-r');

end
