function xt1 = uncontrolled_process(xt0,eps,s,s2,k,h)
%to calculate discretized uncontrolled process x(t+1) = x(t) + fdt + hdw

    xt1 = zeros(3,3);
    xt1(1:2,1:2) = -s2(3)*xt0(1:2,1:2)*h + s(3)*xt0(1:2,1:2)*wedge(eps(3)*sqrt(h));
    xt1(1:2,3) = xt0(1:2,3) + k*xt(1:2,3)*h + s(1)*xt0(1:2,1:2)*eps(1:2)*sqrt(h);

end

