clearvars
close all
run('parameters_SE2.m')

fail_cnt = 0; %number of trajectories failed

for traj_itr = 1:traj_num
    traj_itr
   
    X = zeros(2,(T/h)+1, 'gpuArray'); %to store all positions of this trajectory
    X(1:2,1) = x0(1:2,3); %stack the initial position

    xt = gpuArray(x0); %start the state from the given initial position
    f_xt = gpuArray(f_nominal(xt, s2, k)); %initial f_nominal
    safe_flag_traj = 1;

    for t = t0:h:T-h % this loop is to find u(t) and x(t) at each time step t => x(t+h) = x(t) + f(x(t)).h + G.u(t).h + Sigma*dw
        eps_t_all = randn(n, runs, 'gpuArray'); %GPU array that stores eps(t) at the start of each sample path starting at time t and state xt
        S_tau_all = zeros(1, runs, 'gpuArray'); %GPU array that stores S(tau) of each sample path starting at time t and state xt
      
        parfor i = 1:runs % this loop is to collect S(tau_i) of each sample path starting at time t and state xt

%             eps_t = gather(eps_t_all(:,i)); %standard normal noise at t
%             eps_t = eps_t_all(:,i) %gpuarray

            %parameters at the start of each sample path tau_i
            xt_prime = gpuArray(xt);
            f_xt_prime = gpuArray(f_xt);
%             eps_t_prime = gpuArray(eps_t); %3*1
%             eps_t_runs = randn(n,(T-t)/h, 'gpuArray');
            eps_t_prime = eps_t_all(:,i);
%             xt_runs = randn(n,((T-t)/h)+1, 'gpuArray');
%             xt_runs(1:2,1) = xt_prime(1:2,3);
%             xt_runs(3,1) = xt_prime(1,1)*sin(0.5)+xt_prime(2,1)*cos(0.5); %cos(\theta\2)
            S_tau = 0; %the cost-to-go of the state dependent cost of a sample path
            safe_flag_tau = 1;
%             idx = 1;
               
        
            for t_prime = t:h:T-h % this loop is to compute S(tau_i)
                
%                 S_tau = S_tau + h*b*(xt_prime(1:2,3).')*xt_prime(1:2,3); %add the state dependent running cost
                cos_half = xt_prime(1,1)*sin(0.5)+xt_prime(2,1)*cos(0.5); 
                S_tau = S_tau + h*b*(xt_prime(1:2,3).')*xt_prime(1:2,3) + h*2*acos(cos_half); %後からベクトル計算可能
                xt_prime = uncontrolled_process(xt_prime, eps_t_prime, s, s2, k, h);
%                 xt_runs(1:2,idx+1) = xt_prime(1:2,3);
%                 xt_runs(3,idx+1) = xt_prime(1,1)*sin(0.5)+xt_prime(2,1)*cos(0.5);
            


                if (((xt_prime(1,3)>=xR) && (xt_prime(1,3)<=xS) && (xt_prime(2,3)>=yR) && (xt_prime(2,3)<=yS)) || ((xt_prime(1,3)<=xP) || (xt_prime(1,3)>=xQ) || (xt_prime(2,3)<=yP) || (xt_prime(2,3)>=yQ)))%if yes means t_prime=t_exit
                    S_tau = S_tau + eta; %add the boundary cost to S_tau
                    safe_flag_tau = 0;
                    break; %end this tau 
                end
                
                eps_t_prime = randn(n,1,'gpuArray'); %standard normal noise at new t_prime. Will be used in the next iteration
%                 f_xt_prime = f_nominal(xt_prime, s2, k);
%                 idx = idx+1;
            
%                 f_xt_prime = zeros(3,3);
%                 f_xt_prime(1:2,1:2) = -s2R.*xt_prime(1:2,1:2);
%                 f_xt_prime(1:2,3) = k.*xt_prime(1:2,1:2)*xt_prime(1:2,3); %f_xt_prime at new t_prime. Will be used in the next iteration                   
            end
            if(safe_flag_tau==1) %if tau has not collided 
                  
                cos_half = xt_prime(1,1)*sin(0.5)+xt_prime(2,1)*cos(0.5);
                S_tau = S_tau + d*(xt_prime(1:2,3).')*xt_prime(1:2,3) + d*2*(acos(complex(cos_half)));
                
%                 S_tau = S_tau + h*b*(norm(xt_runs(1,1:end-1))^2 + norm(xt_runs(2,1:end-1))^2) + h*sum(2*acos(complex(xt_runs(3,1:end-1)))) + d*(xt_prime(1:2,3).')*xt_prime(1:2,3) + d*2*(acos(complex(cos_half))); %add the terminal cost to S_tau
%             else
%                 S_tau = S_tau + eta;
%                 S_tau = S_tau + h*b*(norm(xt_runs(1,:))^2 + norm(xt_runs(2,:))^2) + h*sum(2*acos(complex(xt_runs(3,:)))) + eta;
            end
            
            S_tau_all(i) = real(S_tau);
            disp(i);
    
        end
        disp(t);
      
        eps_t_all_arr = gather(eps_t_all); %convert from GPU array to normal array (size: (n X runs))
        S_tau_all_arr = gather(S_tau_all); %convert from GPU array to normal array (size: (1 X runs))

%         denom_i = exp(-S_tau_all_arr/lambda(1)); %(size: (1 X runs))

        denom_iR = exp(-S_tau_all_arr/lambda(3));
        denom_ip = exp(-S_tau_all_arr/lambda(1));

%         numer = eps_t_all_arr*(denom_i.'); %(size: (3 X 1))
        numer_R = wedge(eps_t_all_arr(3,:)*(denom_iR.')/sqrt(h)); %2 X 2
        numer_p = eps_t_all_arr(1:2,:)*(denom_ip.');

%         denom = sum(denom_i); %scalar
        denom_R = sum(denom_iR);
        denom_p = sum(denom_ip);
        
        ut_R = -(s(3)*numer_R)/(denom_R);
        ut_p = -(s(1)*numer_p(1:2,1))/(sqrt(h)*denom_p);

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
        Rt = xt(1:2,1:2);
        xt(1:2,1:2) = Rt + f_xt(1:2,1:2)*h + Rt*ut_R*h + s(3)*Rt*wedge(eps(3)*sqrt(h));
        xt(1:2,3) = xt(1:2,3) + f_xt(1:2,3)*h + Rt*ut_p*h + s(1)*Rt*eps(1:2)*sqrt(h);
%         xt = xt + f_xt*h + xt*(ut_hat)*h + s*xt*eps*sqrt(h); %update the position with the control input ut=> x(t+h) = x(t) + f.h + g.u(t).h + sigma*dw
%         if (xt(1,3) ~= 0) && (xt(2,3) ~= 0)
%         X = [X,xt(1:2,3)]; %stack the new position
        X(1:2,t*100+2) = xt(1:2,3);
%         end
        
        if(((xt(1,3)>=xR) && (xt(1,3)<=xS) && (xt(2,3)>=yR) && (xt(2,3)<=yS)) || ((xt(1,3)<=xP) || (xt(1,3)>=xQ) || (xt(2,3)<=yP) || (xt(2,3)>=yQ))) %if yes means trajectory has crossed the safe set
            fail_cnt = fail_cnt+1; 
            safe_flag_traj = 0;
            break;  %end this traj    
        end 
        
%         f_xt = k*c*xt; %update f(x(t)) for the next t => t=t+h. Will be used in the next iteration
%         f_xt(1:2,1:2) = -s2R.*xt(1:2,1:2);
%         f_xt(1:2,3) = k.*xt(1:2,1:2)*xt(1:2,3);
        f_xt = f_nominal(xt, s2, k);
        X(1:2,1:(t*100)+2)
    end
  
    plot (X(1, 1:(t*100+2)), X(2, 1:(t*100+2)), 'b', 'LineWidth',1)
end
 
fail_prob = fail_cnt/traj_num

figname = ['eta=',num2str(eta),'.fig'];
saveas(gcf,figname)