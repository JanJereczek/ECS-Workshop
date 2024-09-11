hasnan(x) = any(isnan, x)
latexify(x) = L"%$x $\,$"
latexticks(x) = (x, latexify.(x))