function f = phi_approx(x)
    f = 1 ./ (1 + exp(-(0.07056*x.^3+1.5976*x)));
end