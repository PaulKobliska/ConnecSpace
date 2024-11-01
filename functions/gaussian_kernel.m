function r = gaussian_kernel(x,y)

% This function set a gaussian kernel for the rate calculation.
% From CBM, june 2011


r = 0.15915494309190 * exp(-0.5*(x.*x + y.*y));