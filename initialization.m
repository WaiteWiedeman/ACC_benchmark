

function [ population ] = initialization(M, N, Problem)

for i = 1 : M
    for j = 1 : N 
        population.Chromosomes(i).Gene(j) = (Problem.ub(j) - Problem.lb(j)) * rand() + Problem.lb(j);
    end
end

