% Define state derivative for simulation purposes
function [rate] = x_dot_x(t, state, dv, k, mv, kv1, kv2)

disturbance = 0;
measurement_noise = 0.01*sin(2*pi*2*t);

%dv = FIS(gene, state(5)-state(9), state(6)-state(10));

u = -( kv1*(state(5)-state(9)) + dv*(state(7)-state(10)) );

% kalman filter
ym = state(6) + measurement_noise;
[xe] = kalman_filter(ym,u,t,k);
state(5) = xe(:,1);
state(6) = xe(:,2);
state(7) = xe(:,3);
state(8) = xe(:,4);

if u < -1
    u = -1;
elseif u > 1
    u = 1;
end

measurement_error = state(2)+measurement_noise - state(6);

rate = [state(3);...%1
state(4);...%2
-k*state(1)+k*state(2)+u;...%3
k*state(1)-k*state(2)+disturbance;...%4
state(7)+measurement_error;...%5 +M(1)*measurement_error;
state(8)+measurement_error;...%6 +measurement_error;
-state(5)+state(6)+u+measurement_error;...%7+M(3)*measurement_error;
state(5)-state(6)+measurement_error;...%8+M(4)*measurement_error
state(10);...%9
-(kv2/mv)*state(9)+(kv1/mv)*(state(5)-state(9))+(dv/mv)*(state(7)-state(10))];%10

end