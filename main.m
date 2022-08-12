clearvars
close all
run('parameters.m')

fail_cnt = 0; %number of trajectories failed

for traj_itr = 1:traj_num
    traj_itr
   
    X = []; %to store all positions of this trajectory
    X = [X,x0]; %stack the initial position
    
    xt = x0; %start the state from the given initial position
    f_xt = k*c*xt; %initial f_xt (-kpx)
    safe_flag_traj = 1;

    for t = t0:h:T-h % this loop is to find u(t) and x(t) at each time step t => x(t+h) = x(t) + f(x(t)).h + G.u(t).h + Sigma*dw
        eps_t_all = randn(n, runs, 'gpuArray'); %GPU array that stores eps(t) at the start of each sample path starting at time t and state xt
        S_tau_all = zeros(1, runs, 'gpuArray'); %GPU array that stores S(tau) of each sample path starting at time t and state xt
      
        parfor i = 1:runs % this loop is to collect S(tau_i) of each sample path starting at time t and state xt

            eps_t = gather(eps_t_all(:,i)); %standard normal noise at t

            %parameters at the start of each sample path tau_i
            xt_prime = xt;
            f_xt_prime = f_xt;
            eps_t_prime = eps_t;
            S_tau = 0; %the cost-to-go of the state dependent cost of a sample path
            safe_flag_tau = 1;
        
            for t_prime = t:h:T-h % this loop is to compute S(tau_i)
                
                S_tau = S_tau + h*b*(xt_prime.')*xt_prime; %add the state dependent running cost
                
                xt_prime = xt_prime + f_xt_prime*h + s*eps_t_prime*sqrt(h); %move tau ahead
               
                if (((xt_prime(1)>=xR) && (xt_prime(1)<=xS) && (xt_prime(2)>=yR) && (xt_prime(2)<=yS)) || ((xt_prime(1)<=xP) || (xt_prime(1)>=xQ) || (xt_prime(2)<=yP) || (xt_prime(2)>=yQ)))%if yes means t_prime=t_exit
                    S_tau = S_tau + eta; %add the boundary cost to S_tau
                    safe_flag_tau = 0;
                    break; %end this tau 
                end
                
                eps_t_prime = randn(n,1); %standard normal noise at new t_prime. Will be used in the next iteration 
                f_xt_prime = k*c*xt_prime; %f_xt_prime at new t_prime. Will be used in the next iteration                   
            end
          
            if(safe_flag_tau==1) %if tau has not collided 
                S_tau = S_tau + d*(xt_prime.')*xt_prime; %add the terminal cost to S_tau
            end
            
            S_tau_all(i) = S_tau;
        end

        disp(t);
      
        eps_t_all_arr = gather(eps_t_all); %convert from GPU array to normal array (size: (n X runs))
        S_tau_all_arr = gather(S_tau_all); %convert from GPU array to normal array (size: (1 X runs))

        denom_i = exp(-S_tau_all_arr/lambda); %(size: (1 X runs))
        numer = eps_t_all_arr*(denom_i.'); %(size: (2 X 1))
        denom = sum(denom_i); %scalar

        ut = ((s/c)*numer)/(sqrt(h)*denom); %the control input
        
        if(any(isnan(ut(:))))
            fprintf("error!")
            return
        end
        
        %move the trajectory forward
        eps = randn(n,1);
        xt = xt + f_xt*h + c*ut*h + s*eps*sqrt(h); %update the position with the control input ut=> x(t+h) = x(t) + f.h + g.u(t).h + sigma*dw
        X = [X,xt]; %stack the new position
        
        if(((xt(1)>=xR) && (xt(1)<=xS) && (xt(2)>=yR) && (xt(2)<=yS)) || ((xt(1)<=xP) || (xt(1)>=xQ) || (xt(2)<=yP) || (xt(2)>=yQ))) %if yes means trajectory has crossed the safe set
            fail_cnt = fail_cnt+1; 
            safe_flag_traj = 0;
            break;  %end this traj    
        end 
        
        f_xt = k*c*xt; %update f(x(t)) for the next t => t=t+h. Will be used in the next iteration 
    end
  
    plot (X(1, :), X(2, :), 'b', 'LineWidth',1)
end
 
fail_prob = fail_cnt/traj_num

figname = ['eta=',num2str(eta),'.fig'];
saveas(gcf,figname)