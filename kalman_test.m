clc; clear; close all;
endtime = 40;
% nominal values
k_n = 1;
m1_n = 1;
m2_n = 1;
mv = 600;
kv1 = 0.35;
kv2 = 2.9;
k_range = [0.5, 2];
k = k_range(1) + (k_range(2)-k_range(1))*rand(1);
% Simulate
x_0 = [0; 0; 0; 1; 0; 0; 0; 0; 0; 0];
u = 1;
%options = odeset('AbsTol',1e-3);
%ode45?
[ts,xs] = ode45(@(t,x) x_dot_x(t, x, u, k, mv, kv1, kv2), [0 endtime], x_0);%, options);
ym = xs(6,:);
m = length(ts);
% Discrete-Time Process Noise Covriance
q=1e-5*eye(10);

% Initial Covariance
poa=1e-5;
p=poa*diag([1 1 1 1 1 1 1 1 1 1]);
pcov=zeros(m,2);pcov(1,:)=poa*[1 1 1 1 1 1 1 1 1 1];

% Initial Condition and H Matrix (constant)
x0=[ym(1);0; 0;0;0;0;0;0;0;0 ];xe=zeros(m,10);xe(1,:)=x0';x=x0; %
h=[0 1 0 0];

% Main Loop
for i = 1:m-1

% Kalman Gain    
gain=p*h'*inv(h*p*h'+poa^2);
%disp(gain)
% Update
x=x+gain*(ym(i)-h*x);
p=[eye(2)-gain*h]*p;
%disp(p)
% Propagate
x=A*x+B*u(m);
p=A*p*A'+q;
%disp(p)
% Store Variables
xe(i+1,:)=x';
pcov(i+1,:)=diag(p)';

end

% 3-Sigma Outlier
sig3=pcov.^(0.5)*3;

figure; 
plot(t,y,'--',t,ym )
figure; 
plot(t,xe(:,1),'--',t,x_true(:,1) )
figure; 
plot(t,xe(:,2),'--',t,x_true(:,2) )
% figure; 
% plot(t,xe(:,3),'--',t,x_true(:,3) )
% figure; 
% plot(t,xe(:,4),'--',t,x_true(:,4) )
% Plot Results
% plot(t,[sig3(:,1) xe(:,1)-x_true(:,1) -sig3(:,1)])
% grid;
% set(gca,'Fontsize',12);set(gca,'Fontname','Helvetica');
% xlabel('Time (Min)');
% hh=get(gca,'Ylabel');
% set(hh,'String','\fontsize{12} {Attitude Error ({\mu}rad)}');
% 
% disp(' Press any key to continue')
% pause
% 
% plot(t,xe(:,2));grid
% set(gca,'Fontsize',12);set(gca,'Fontname','Helvetica');
% xlabel('Time (Min)');
% ylabel('x2 est')