function [p_hat] = upper_bound_confidence_binomial_approximation(c, n, confidence)
% Approximate upper confidence bound

epsilon = 1 - confidence;

% approximate upper bound of poisson mean
%poisson_cum_error(8.4060, 2, 1-.99)
% options = optimoptions('fsolve');
% options = optimoptions(options, 'FunctionTolerance', (5e-7)^4, 'OptimalityTolerance', 1e-21);
% [m_hat, err] = fsolve(@(m) poisson_cum_error(m, c, epsilon), c, options);

m_hat = 0.5 * chi2inv(epsilon, 2*c);

% determine value of f
if confidence == 0.99
    if (0 <= c) && (c <= 2)
        f = 0.55;
    elseif (3 <= c) && (c <= 33)
        f = .75;
    else
        f = 0;
    end
elseif confidence == 0.95
    if (0 <= c) && (c <= 21)
        f = 0.45;
    else
        f = 0;
    end
else
    error('Invalid confidence value')
end

% approximate upper bound of binomial USING THE 2ND APPROXIMATION
p_hat = (m_hat/n) * ( n - 0.25*m_hat + f ) / ( n + 0.25*(m_hat - 2*c) + f );
end


function [err] = poisson_cum_error(m, r, epsilon)
err = 0;
for x = 0:r
    err = err + ((m^x)/factorial(x));
end

err = err * exp(-m);
err = err - epsilon;
end