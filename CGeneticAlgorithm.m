
function [BestChrom, cgcurve]  = CGeneticAlgorithm (M , N, MaxGen , Pc, Pm , Er , Problem, order)

cgcurve = zeros(1,MaxGen);

%%  Initialization
g = 1;
disp(['Generation #' , num2str(g)]);
[ population ] = initialization(M, N, Problem);
if order
    population = ordering(population);
end


%%%% CALCULATE FITNESS
% disp('duplicate')
% tic
% somes = population.Chromosomes;
% objectivefunc = @(x) Problem.obj(x);
% toc
% 
% disp('calc fitness')
% tic
% for i = 1 : M
%     somes(i).fitness = objectivefunc( somes(i).Gene(:) );
% end
% toc
% 
% disp('replace')
% tic
% population.Chromosomes = somes;
% toc
population.Chromosomes = calc_pop_fit(population, Problem, M);
%%%%

all_fitness_values = [ population.Chromosomes(:).fitness ];
[cgcurve(1) , ~ ] = max( all_fitness_values);
disp(['Max Fitness' , num2str(cgcurve(1))]);



%% Main loop
for g = 2:MaxGen
    
    disp(['Generation #' , num2str(g)]);
%     % Calcualte the fitness values
%     for i = 1 : M
%         population.Chromosomes(i).fitness = Problem.obj( population.Chromosomes(i).Gene(:) );
%     end
    
    %drawnow
    
    for k = 1: 2: M
        % Selection
        [ parent1, parent2] = selection(population);
        
        % Crossover
        [child1 , child2] = crossover_continious(parent1 , parent2, Pc, Problem);
        
        % Mutation
        
        %[child1] = mutation_continious(child1, Pm, Problem);
        %[child2] = mutation_continious(child2, Pm, Problem);

        [child1] = mutation_continious_local(child1, Pm, Problem);
        [child2] = mutation_continious_local(child2, Pm, Problem);
        
        newPopulation.Chromosomes(k).Gene = child1.Gene;
        newPopulation.Chromosomes(k+1).Gene = child2.Gene;
    end
    
    if order
    	newPopulation = ordering(newPopulation);
    end
    
    % Elitism
    [ newPopulation ] = elitism(population , newPopulation, Er);

    population = newPopulation;

    % determine fitness of population
%     for i = 1 : M
%         population.Chromosomes(i).fitness = Problem.obj( population.Chromosomes(i).Gene(:) );
%     end
    population.Chromosomes = calc_pop_fit(population, Problem, M);

    
    all_fitness_values = [ population.Chromosomes(:).fitness ];
    [cgcurve(g) , ~ ] = max(all_fitness_values);
    %cgcurve(g)
    disp(['Max Fitness: ' , num2str(cgcurve(g))]);
end

% for i = 1 : M
%     population.Chromosomes(i).fitness = Problem.obj( population.Chromosomes(i).Gene(:) );
% end


[max_val , indx] = sort([ population.Chromosomes(:).fitness ] , 'descend');
    
BestChrom.Gene    = population.Chromosomes(indx(1)).Gene;
BestChrom.Fitness = population.Chromosomes(indx(1)).fitness;

% figure(1)
% plot(cgcurve)
% hold on
 
end