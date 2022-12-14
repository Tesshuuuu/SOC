function f_xt = f_nominal(xt, s2, k)
%to calculate function f(t)

    f_xt = zeros(3,3);
    f_xt(1:2,1:2) = -s2(3)*xt(1:2,1:2);
    f_xt(1:2,3) = k*xt(1:2,3);
end

