
function [out]  = FIS(gene, x, xd)
gene(1:5) = sort(gene(1:5));
gene(6:10) = sort(gene(6:10));
% extract centers of MF's from gene and evaluate MF's
x_c1 = gene(1);
x_c2 = gene(2);
x_c3 = gene(3);
x_c4 = gene(4);
x_c5 = gene(5);

xd_c1 = gene(6);
xd_c2 = gene(7);
xd_c3 = gene(8);
xd_c4 = gene(9);
xd_c5 = gene(10);

% determine membership values
mu_x1 = lshlder(x, x_c1, x_c2, -2);
mu_x2 = triangle(x, x_c2, x_c3, x_c1);
mu_x3 = triangle(x, x_c3, x_c4, x_c2);
mu_x4 = triangle(x, x_c4, x_c5, x_c3);
mu_x5 = rshlder(x, x_c5, 2, x_c4);

mu_xd1 = lshlder(xd, xd_c1, xd_c2, -2);
mu_xd2 = triangle(xd, xd_c2, xd_c3, xd_c1);
mu_xd3 = triangle(xd, xd_c3, xd_c4, xd_c2);
mu_xd4 = triangle(xd, xd_c4, xd_c5, xd_c3);
mu_xd5 = rshlder(xd, xd_c5, 2, xd_c4);

% Find antecedent memberships using multiplication of the input memberships
rule_memberships = [mu_xd1; mu_xd2; mu_xd3; mu_xd4; mu_xd5] * [mu_x1, mu_x2, mu_x3, mu_x4, mu_x5];
rule_memberships = reshape(rule_memberships,1,25);

% Evaluate consequents of TKS
rule_vals = gene(11:35);
rule_vals = reshape(rule_vals,1,25);
%disp(sum(rule_memberships))
%disp('rule_memberships')
%disp(rule_memberships)
%disp('rule_vals')
%disp(rule_vals)

out = sum(rule_memberships .* rule_vals);%/sum(rule_memberships);
%disp('out')
%disp(out)
end