function [p_hat] = lower_bound_confidence_binomial_approximation(c, n, confidence)
% Approximate upper confidence bound

epsilon = 1 - confidence;

% approximate upper bound of poisson mean
% poisson_cum_error(5.4282, 12, 1-.99)
% options = optimoptions('fsolve');
% options = optimoptions(options, 'FunctionTolerance', (5e-7)^4, 'OptimalityTolerance', 1e-21);
% [m_hat, err] = fsolve(@(m) poisson_cum_error(m, c, epsilon), c, options)

m_hat = 0.5 * chi2inv(epsilon, 2*c);

% determine value of g
if confidence == 0.99
    if (10 <= c) && (c <= 45)
        g = 0.70;
    elseif c <= 9
        g = 0.070 * c;
    else
        g = 0;
    end
elseif confidence == 0.95
    if (10 <= c) && (c <= 28)
        g = 0.60;
    elseif c <= 9
        g = 0.060 * c;
    else
        g = 0;
    end
else
    error('Invalid confidence value')
end

% approximate upper bound of binomial USING THE 2ND APPROXIMATION
p_hat = (m_hat/n) * ( n - m_hat/3 - g ) / ( n - 0.5*(c-1) + m_hat/6 - g );
end


function [err] = poisson_cum_error(m, r, epsilon)
err = 0;
for x = 0:r
    err = err + ((m^x)/factorial(x));
end

err = err * exp(-m);
err = err - epsilon;
end