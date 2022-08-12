function [xt] = rollout(xt, eps)
%to compute sequence of xt(size:n*T-h)

    t=0;
%     seq = size(xt,2);
%     xt_prime = [cos(xt(3)) -sin(xt(3)) xt(1); sin(xt(3)) cos(xt(3)) xt(2); 0 0 1];
    while t<=T-t
%         S_tau = state_cost(S_tau, h, b, xt_prime); %add the state dependent running cost
        t = t+h;
        xt_prime = uncontrolled_process(xt_prime, eps(:,t/h+1));
        xt(1:2,t/h+1) = xt_prime(1:2,3);
        xt(3,t/h+1) = calculate_theta(xt_prime);
        if (((xt_prime(1,3)>=xR) && (xt_prime(1,3)<=xS) && (xt_prime(2,3)>=yR) && (xt_prime(2,3)<=yS)) || ((xt_prime(1,3)<=xP) || (xt_prime(1,3)>=xQ) || (xt_prime(2,3)<=yP) || (xt_prime(2,3)>=yQ)))%if yes means t_prime=t_exit
%             S_tau = S_tau + eta; %add the boundary cost to S_tau
%             safe_flag_tau = 0;
            break; %end this tau 
        end
%         eps_t_prime = randn(n,1);
%         f_xt_prime = f_nominal(xt_prime, s2, k);
    end

end

function S1 = state_cost(S0, h, b, x0)
    S1 = S0 + h*b*(x0(1:2,3).')*x0(1:2,3);
end

function x1 = uncontrolled_process(x0, f, eps)
    x1(1) = -s2(3)*xt(1)*h + s(3)*x0(1)*

    x1(1:2,1:2) = f(1:2,1:2)*h + s(3)*x0(1:2,1:2)*wedge(eps(3)*sqrt(h));
    x1(1:2,3) = f(1:2,3)*h + s(1:2)'*x0(1:2,1:2)*eps(1:2)*sqrt(h);
end

function f_xt = f_nominal(xt)
%to calculate uncontrolled process x(t+1) = x(t) + fdt + hdw

    f_xt(1) = -s2(3)
    f_xt(1:2,1:2) = -s2(3)*xt(1:2,1:2);
    f_xt(1:2,3) = k*xt(1:2,1:2)*xt(1:2,3);

end

function theta = calculate_theta(x)
%CALCULATE_THETA return theta from rotation matrix
    if size(x,1) == 3
        theta_cand = zeros(2,1);
        theta_cand(1) = asin(x(2,1));
        theta_cand(2) = sign(theta_cand(1))*pi - theta_cand(1);
        if fix(cos(theta_cand(1))*10) == fix(x(1,1)*10)
            theta = theta_cand(1)
        else
            theta = theta_cand(2)
        end
    end

end