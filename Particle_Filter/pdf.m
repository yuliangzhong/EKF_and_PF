function pw = pdf(w,eps)
    w = abs(w);
    if w < 2*eps
        pw = -w/(5*eps^2) + 2/5/eps;
    elseif w<2.5*eps
        pw = 2/(5*eps^2)*(w - 2*eps);
    elseif w<3*eps
        pw = -2/(5*eps^2)*(w - 3*eps);
    else
        pw = 0;
    end
end