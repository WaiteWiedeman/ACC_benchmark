
%****************************************************************************
 % A 'GENERIC' ALGORITHM IN MATLAB CODE THAT RETURNS DEFUZZIFIED Decision Applied
 
 % to MATLAB's Tipper Problem. THIS PROGRAM WAS WRITTEN BY KELLY COHEN WITHIN THE FRAMEWORK 
 
 % The SOFT COMPUTING BASED AI Class Taught at UC/CEAS
 %               -------------------------------------
 %             -  LAST UPDATE : September 14, 2019     -
 %               -------------------------------------
 %****************************************************************************
%
function mu = triangle(x,ce,re,le)    %TRIANGLE MEMBERSHIP FUNCTION
if x >= le && x < ce
    if le == ce
    mu=1;
    else
    mu = (x - le)/(ce - le);
    end
elseif x >= ce && x <= re
    if re == ce
    mu=1;
    else
    mu = (re - x)/(re - ce);
    end
else
mu = 0;
end

