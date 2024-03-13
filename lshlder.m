%****************************************************************************
 % A 'GENERIC' ALGORITHM IN MATLAB CODE THAT RETURNS DEFUZZIFIED Decision Applied
 
 % to MATLAB's Tipper Problem. THIS PROGRAM WAS WRITTEN BY KELLY COHEN WITHIN THE FRAMEWORK 
 
 % The SOFT COMPUTING BASED AI Class Taught at UC/CEAS
 %               -------------------------------------
 %             -  LAST UPDATE : September 14, 2019     -
 %               -------------------------------------
 %****************************************************************************
%
function mul = lshlder(x,cel,rel,lel)   %left shoulder membership function


if x>cel && x<=rel
    if rel == cel
    mul=1;
    else
    mul = (rel-x)/(rel-cel);
    end
elseif x<=cel
mul = 1;
else
mul = 0;
end

