function y = gaussian(x,sig,c)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    y = exp(- (x-c).^2 ./ (2*sig.^2) );
end

