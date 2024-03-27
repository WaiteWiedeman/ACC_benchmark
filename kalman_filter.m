function [xe] = kalman_filter(ym,u,t,k)
m=length(t);
m1_n = 1;
m2_n = 1;
A=[0 0 1 0; 0 0 0 1; -k/m1_n k/m1_n 0 0; k/m2_n -k/m2_n 0 0];
B=[1 0 0 0]';
C=[0 1 0 0];
% Discrete-Time Process Noise Covriance
q=eye(4);

% Initial Covariance
poa=1e-4;
p=poa*eye(4);
pcov=zeros(m,4);pcov(1,:)=[poa poa poa poa];

% Initial Condition and H Matrix (constant)
x0=[0;ym(1);0;0];xe=zeros(m,4);xe(1,:)=x0';x=x0;
h=C;

% Main Loop
for i = 1:m-1

% Kalman Gain    
gain=p*h'*inv(h*p*h'+poa^2);

% Update
x=x+gain*(ym(i)-h*x);
p=[eye(4)-gain*h]*p;

% Propagate
x=A*x+B*u;
p=A*p*A'+q;

% Store Variables
xe(i+1,:)=x';
pcov(i+1,:)=diag(p)';

end
