function theta = calculate_theta(x)
%CALCULATE_THETA return theta from rotation matrix
    if size(x,1) == 3
        theta_cand = zeros(2,1);
        theta_cand(1) = asin(x(2,1));
        theta_cand(2) = sign(theta_cand(1))*pi - theta_cand(1);
        if fix(cos(theta_cand(1))*10) == fix(x(1,1)*10)
            theta = theta_cand(1);
        else
            theta = theta_cand(2);
        end
    end

end

