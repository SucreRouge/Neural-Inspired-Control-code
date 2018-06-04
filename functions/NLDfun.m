function y = NLDfun(s,NLDg,NLDs)
%#codegen

% u = NLDfun(s)
% NLDs = [0.5];
% NLDg = [10];
posNeg = sign(s);
y = (  ( 1./ (1+ exp(-NLDg.*(abs(s)-NLDs)) ) -0.5) +0.5).*posNeg; 