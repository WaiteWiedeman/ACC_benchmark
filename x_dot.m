% Define state derivative for simulation purposes
function [rate] = x_dot(t, state, gene, M, k, mv, kv1, kv2)

disturbance = 0;
measurement_noise = 0; %0.01*sin(2*pi*2*t);

%dv = FIS(gene, state(5)-state(9), state(6)-state(10));
FIS1 = FIS(gene(1:35),state(5),state(3));
FIS2 = FIS(gene(36:70),state(6),state(4));
dv = FIS(gene(71:105), FIS1, FIS2);

u = -( kv1*(state(5)-state(9)) + dv*(state(7)-state(10)) );

% disp('dv')
% disp(dv)
% disp('u')
% disp(u)
% disp('state')
% disp(state)
% disp('gene')
% disp(gene)
% disp('state5')
% disp(state(5))
% disp('state6')
% disp(state(6))

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
state(7)+M(1)*measurement_error;...%5
state(8)+M(2)*measurement_error;...%6
-state(5)+state(6)+u+M(3)*measurement_error;...%7
state(5)-state(6)+M(4)*measurement_error;...%8
state(10);...%9
-(kv2/mv)*state(9)+(kv1/mv)*(state(5)-state(9))+(dv/mv)*(state(7)-state(10))];%10

end