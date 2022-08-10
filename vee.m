function vec = vee(mat)
    if size(Rt,1) == 2
        vec = mat(2,1);
    elseif size(Rt,1) == 3
        vec = zeros(3,1);
        vec(1) = mat(3,2);
        vec(2) = mat(1,3);
        vec(3) = mat(2,1);
    end
end