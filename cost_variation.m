%% Plot variation of cost with q
qs = linspace(120,170,100);

costs = zeros(size(qs));

for m = 1:length(costs)
q = qs(m) % psf

Ks = FIS(BestChrom.Gene, q);
K1 = Ks(1);
K2 = Ks(2);
K3 = Ks(3);
F = [0 1; -K2 -K1];
G = diag([K1 K2])*[0 1 0 0 0 0; 0 0 0 1 0 0];
F_p = [G F];

% Define BACT model for given q, U0
[Ab, Bb, AaTE, BaTE, CaTE, DaTE, Ag, Bg, Cg, Dg, Ep] = bact_model(q, U0, 0);
A = [Ab Bb*CaTE; zeros(2,4) AaTE];
B = [Bb*DaTE; BaTE];

% Simulate
z_0 = [0; 0; 0; 10*pi/180; 0; 0; 0; 0];
[ts,xs] = ode45(@(t,x)xdot_passive(A,B,F_p,K3,x,t), [0 5], z_0);

ubound = 12*pi/180;
cost = 0;
for n = 1:length(ts)
    state = xs(n,:);
    yd = F_p*state';
    u(n) = min(max(-K3*yd(2),-ubound),ubound);
    un = u(n);
    
    cost = cost + (state(3))^2 + (state(4))^2 + 2*(un^2);
end
cost = cost/length(ts);
costs(m) = cost;

end

% Plot results
figure(1)
plot(qs,costs)
xlabel('q (psf)')
ylabel('Cost')













