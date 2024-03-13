
function [population] = ordering(population)

M = length(population.Chromosomes); %length(population.Chromosomes(:));

% just for ensuring that membership functions centers are sorted after
% initialization and mutation
for i = 1 : M
    x_c = population.Chromosomes(i).Gene(1:5);
    y_c = population.Chromosomes(i).Gene(6:10);
    
    population.Chromosomes(i).Gene(1:5) = sort(x_c);
    population.Chromosomes(i).Gene(6:10) = sort(y_c);
end

end