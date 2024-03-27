
function [fit]  = acc_fitness(gene, N_monte_carlo, endtime, M, k_range, mv, kv1, kv2)
% seed rng
rng(1)

u_above_thresh = zeros(1,N_monte_carlo);
settling_time_above_thresh = zeros(1,N_monte_carlo);
x_above_thresh = zeros(1,N_monte_carlo);
%u2tot = 0;
tscnt = 0;
for sim_n = 1:N_monte_carlo
    u = [];
    % determine k value
    k = k_range(1) + (k_range(2)-k_range(1))*rand(1);

    % Simulate
    x_0 = [0; 0; 0; 1; 0; 0; 0; 0; 0; 0];
    %options = odeset('AbsTol',1e-3);
    %ode45?
    [ts,xs] = ode23(@(t,x) x_dot(t, x, gene, M, k, mv, kv1, kv2), [0 endtime], x_0);%, options); 
    info1 = lsiminfo(xs(:,1),ts,'SettlingTimeThreshold',0.1);
    T_s1(sim_n) = info1.TransientTime;
    info2 = lsiminfo(xs(:,2),ts,'SettlingTimeThreshold',0.1);
    T_s2(sim_n) = info2.TransientTime;
    
    % Calculate cost
    for n = 1:length(ts)
        tscnt = tscnt + 1;
        t = ts(n);
        state = xs(n,:);

        % u exceeds bounds
        %dv = FIS(gene, state(5)-state(9), state(6)-state(10));
        FIS1 = FIS(gene(1:35),state(5),state(7));
        FIS2 = FIS(gene(36:70),state(6),state(8));
        dv = FIS(gene(71:105), FIS1, FIS2);

        u(n) = -( kv1*(state(5)-state(9)) + dv*(state(7)-state(10)) );
        %u2tot = u2tot + u^2;
        if abs(u(n)) > 1
            u_above_thresh(sim_n) = 1;
        end

        max_state = max(abs(state(1)), abs(state(2)));
        % settling time > 15 seconds
        if t >= 15
            if max_state > 0.1
                settling_time_above_thresh(sim_n) = 1;
            end
        end

        % x above bounds ('unstable')
        if max_state > 2
            x_above_thresh(sim_n) = 1;
        end     

        %cost = cost + (state(3))^2 + (state(4))^2 + 2*(u^2);
    end
    us(sim_n) = max(u);
end

% disp(sum(x_above_thresh)/N_monte_carlo)
% disp(sum(settling_time_above_thresh)/N_monte_carlo)
% disp(sum(u_above_thresh)/N_monte_carlo)
% disp(sum(us))
% disp(mean(T_s1)+mean(T_s2))

% fitness is weighted squared probability
fit = (1e5)*(1*(sum(x_above_thresh)/N_monte_carlo)^2 + 0.01*(sum(settling_time_above_thresh)/N_monte_carlo)^2 ... 
    + 0.1*(sum(u_above_thresh)/N_monte_carlo)^2 + 0.001*sum(us)^2 + 0.00001*(mean(T_s1)+mean(T_s2))^2);
% negative because fitness
%fit = -fit;

% prefer less control authority, improve convergence
%fit = fit - (1e0 * u2tot/tscnt);
end
