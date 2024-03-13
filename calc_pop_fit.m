function [somes] = calc_pop_fit(population, Problem, M)
%disp('duplicate')
%tic
somes = population.Chromosomes;
objectivefunc = @(x) Problem.obj(x);
%toc

disp('Calculating fitness')
tic
parfor i = 1 : M
    somes(i).fitness = objectivefunc( somes(i).Gene(:) );
end
toc

end

