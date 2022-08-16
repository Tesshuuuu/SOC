function M = calc_M(eps, sigma, k, h)
%CALC_M calculate 3x3 matix M
    M = zeros(3,3);
    M(1:2,1:2) = (1 - sigma^2)*eye(2) + sigma*wedge(eps(3)*sqrt(h));
    M(1:2,3) = sigma*sqrt(h)*eps(1:2);
    M(3,1:2) = zeros(1,2);
    M(3,3) = 1 - k*h;
end

