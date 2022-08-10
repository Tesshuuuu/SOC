function mat = wedge(vec)
    if length(vec) == 1
        mat = zeros(2,2);
        mat(1,2) = -vec;
        mat(2,1) = vec;
    elseif length(vec) == 3
        mat = [0 -vec(3) vec(2); vec(3) 0 -vec(1); -vec(2) vec(1) 0];
    end
end