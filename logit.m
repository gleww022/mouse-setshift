function y = logit(x)
    y = exp(x)./(1+exp(x));
end