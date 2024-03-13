

function [child] = mutation_continious_local(child, Pm, Problem)

Gene_no = length(child.Gene);
mag_scale = 1/3;

for k = 1: Gene_no
    R = rand();
    mag = mag_scale * (Problem.ub(k) - Problem.lb(k));
    if R < Pm
        child.Gene(k) = child.Gene(k) + mag * randn(1); % (Problem.ub(k) - Problem.lb(k)) * rand() + Problem.lb(k);
        child.Gene(k) = max(child.Gene(k), Problem.lb(k));
        child.Gene(k) = min(child.Gene(k), Problem.ub(k));
    end
end

end