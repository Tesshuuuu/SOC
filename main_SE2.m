clearvars
close all
run('parameters_SE2.m')

fail_cnt = 0; %number of trajectories failed

for traj_itr = 1:traj_num
    traj_itr
   
    X = []; %to store all positions of this trajectory
    X = [X,x0(1:2,3)]; %stack the initial position
    
    xt = x0; %start the state from the given initial position
%     f_xt = k*c*xt; %initial f_xt (-kpx)
    f_xt = zeros(3,3);
    f_xt(1:2,1:2) = -s2R.*xt(1:2,1:2);
    f_xt(1:2,3) = k.*xt(1:2,1:2)*xt(1:2,3);
    safe_flag_traj = 1;

    for t = t0:h:T-h % this loop is to find u(t) and x(t) at each time step t => x(t+h) = x(t) + f(x(t)).h + G.u(t).h + Sigma*dw
        eps_t_all = randn(n, runs, 'gpuArray'); %GPU array that stores eps(t) at the start of each sample path starting at time t and state xt
        S_tau_all = zeros(1, runs, 'gpuArray'); %GPU array that stores S(tau) of each sample path starting at time t and state xt
      
        parfor i = 1:runs % this loop is to collect S(tau_i) of each sample path starting at time t and state xt

            eps_t = gather(eps_t_all(:,i)); %standard normal noise at t

            %parameters at the start of each sample path tau_i
            xt_prime = xt;
            f_xt_prime = f_xt;
            eps_t_prime = eps_t; %3*1
            S_tau = 0; %the cost-to-go of the state dependent cost of a sample path
            safe_flag_tau = 1;
        
            for t_prime = t:h:T-h % this loop is to compute S(tau_i)
                
                S_tau = S_tau + h*b*(xt_prime(1:2,3).')*xt_prime(1:2,3); %add the state dependent running cost
                
                Rt_prime = xt_prime(1:2,1:2);
                xt_prime = xt_prime + f_xt_prime*h; %move tau ahead
                xt_prime(1:2,1:2) = xt_prime(1:2,1:2) + sR*Rt_prime*wedge(eps_t_prime(3)*sqrt(h));
                xt_prime(1:2,3) = xt_prime(1:2,3) + s*Rt_prime*eps_t_prime(1:2)*sqrt(h);

                if (((xt_prime(1,3)>=xR) && (xt_prime(1,3)<=xS) && (xt_prime(2,3)>=yR) && (xt_prime(2,3)<=yS)) || ((xt_prime(1,3)<=xP) || (xt_prime(1,3)>=xQ) || (xt_prime(2,3)<=yP) || (xt_prime(2,3)>=yQ)))%if yes means t_prime=t_exit
                    S_tau = S_tau + eta; %add the boundary cost to S_tau
                    safe_flag_tau = 0;
                    break; %end this tau 
                end
                
                eps_t_prime = randn(n,1); %standard normal noise at new t_prime. Will be used in the next iteration 
                f_xt_prime = zeros(3,3);
                f_xt_prime(1:2,1:2) = -s2R.*xt_prime(1:2,1:2);
                f_xt_prime(1:2,3) = k.*xt_prime(1:2,1:2)*xt_prime(1:2,3); %f_xt_prime at new t_prime. Will be used in the next iteration                   
            end
          
            if(safe_flag_tau==1) %if tau has not collided 
                S_tau = S_tau + d*(xt_prime(1:2,3).')*xt_prime(1:2,3); %add the terminal cost to S_tau
            end
            
            S_tau_all(i) = S_tau;
        end
      
        eps_t_all_arr = gather(eps_t_all); %convert from GPU array to normal array (size: (n X runs))
        S_tau_all_arr = gather(S_tau_all); %convert from GPU array to normal array (size: (1 X runs))

        denom_i = exp(-S_tau_all_arr/lambda); %(size: (1 X runs))
        
        numer = eps_t_all_arr*(denom_i.'); %(size: (3 X 1))
        denom = sum(denom_i); %scalar
        
        ut_R = -(sR*numer(3,1))/(sqrt(h)*denom);
        ut_p = -(s*numer(1:2,1))/(sqrt(h)*denom);

%         ut = -(s*numer)/(sqrt(h)*denom); %the control input
        
        if(any(isnan(ut_p)))
            fprintf("error!")
            return
        end
        if(any(isnan(ut_R)))
            fprintf("error!")
            return
        end



        %move the trajectory forward
        eps = randn(n,1);
%         ut_hat = zeros(3,3);
%         ut_hat(1:2,1:2) = wedge(ut(3));
%         ut_hat(1:2,3) = ut(1:2);
        xt(1:2,1:2) = xt(1:2,1:2) + f_xt(1:2,1:2)*h + xt(1:2,1:2)*wedge(ut_R)*h + sR*xt(1:2,1:2)*wedge(eps(3)*sqrt(h));
        xt(1:2,3) = xt(1:2,3) + f_xt(1:2,3)*h + xt(1:2,1:2)*ut_p*h + s*xt(1:2,1:2)*eps(1:2)*sqrt(h);
%         xt = xt + f_xt*h + xt*(ut_hat)*h + s*xt*eps*sqrt(h); %update the position with the control input ut=> x(t+h) = x(t) + f.h + g.u(t).h + sigma*dw
        X = [X,xt(1:2,3)]; %stack the new position
        
        if(((xt(1,3)>=xR) && (xt(1,3)<=xS) && (xt(2,3)>=yR) && (xt(2,3)<=yS)) || ((xt(1,3)<=xP) || (xt(1,3)>=xQ) || (xt(2,3)<=yP) || (xt(2,3)>=yQ))) %if yes means trajectory has crossed the safe set
            fail_cnt = fail_cnt+1; 
            safe_flag_traj = 0;
            break;  %end this traj    
        end 
        
%         f_xt = k*c*xt; %update f(x(t)) for the next t => t=t+h. Will be used in the next iteration
        f_xt(1:2,1:2) = -s2R.*xt(1:2,1:2);
        f_xt(1:2,3) = k.*xt(1:2,1:2)*xt(1:2,3);
    end
  
    plot (X(1, :), X(2, :), 'b', 'LineWidth',1)
end
 
fail_prob = fail_cnt/traj_num

figname = ['eta=',num2str(eta),'.fig'];
saveas(gcf,figname)