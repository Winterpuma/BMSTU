function lab2()
    X = [3.38,1.21,1.85,2.24,4.17,2.99,4.81,2.71,2.70,4.41,3.21,3.15,2.77,4.05,3.89,1.56,2.78,2.04,2.82,3.28,2.63,1.89,3.57,3.15,3.80,5.40,3.25,2.04,2.61,5.06,2.87,2.66,4.80,3.86,0.09,2.45,2.40,2.14,1.69,2.36,5.44,2.77,1.94,2.55,3.97,1.88,3.01,4.21,4.74,2.02,2.38,2.46,3.51,2.89,1.57,3.53,0.77,3.31,3.58,2.77,3.61,3.71,2.38,3.06,4.29,4.76,1.69,1.59,3.21,2.74,3.99,3.53,3.52,2.84,1.21,2.82,4.34,3.65,2.22,2.87,3.14,3.58,1.96,3.41,3.85,1.96,3.02,4.22,3.10,2.68,3.67,1.70,5.47,5.02,2.52,3.09,2.19,4.44,2.33,2.27,3.34,3.05,4.35,3.58,3.43,4.49,3.57,3.20,1.53,3.53,3.53,1.27,3.40,4.53,2.21,3.28,3.50,2.01,3.30,1.86];
      
    N = length(X);
    gamma = 0.9;
    
    alpha = (1 - gamma) / 2;
    
    mu = get_mu(X);
    Ssqr = get_Ssqr(X);
    array_of_mu = get_array_of_mu(X, N);
    array_of_var = get_array_of_var(X, N);
    mu_higher = get_mu_higher(array_of_mu, array_of_var, alpha, N);
    mu_lower = get_mu_lower(array_of_mu, array_of_var, alpha, N);
    sigma_higher = get_sigma_higher(array_of_var, alpha, N);
    sigma_lower = get_sigma_lower(array_of_var, alpha, N);
    
    
    fprintf('mu = %.3f\n', mu);
    fprintf('S^2 = %.3f\n', Ssqr);
    fprintf('mu_lower = %.3f\n', mu_lower(end));
    fprintf('mu_higher = %.3f\n', mu_higher(end));
    fprintf('sigma^2_lower = %.3f\n', sigma_lower(end));
    fprintf('sigma^2_higher = %.3f\n', sigma_higher(end));
    

    figure 
    hold on;
    plot([1, N], [mu, mu], 'g');
    plot((1 : N), array_of_mu, 'r');
    plot((1 : N), mu_lower, 'b');
    plot((1 : N), mu_higher, 'm');
    legend('\mu\^(x_N)','\mu\^(x_n)','_{--}\mu^(x_n)','^{--}\mu^(x_n)');
    grid on;
    hold off;
    
    figure 
    hold on;
    plot([1, N], [Ssqr, Ssqr], 'g');
    plot((1 : N), array_of_var, 'r');
    plot((1 : N), sigma_lower, 'b');
    plot((4 : N), sigma_higher(4 : length(sigma_higher)), 'm');
    legend('S^2(x_N)','S^2(x_n)','_{--}\sigma^2(x_n)','^{--}\sigma^2(x_n)');
    grid on;
    hold off;
    
    
    function mu = get_mu(X)
        mu = sum(X) / size(X, 2);
    end
    
    function Ssqr = get_Ssqr(X)
        n = size(X, 2);
        Ssqr = n / (n - 1) * get_sigma_sqr(X);
    end

    function array_of_mu = get_array_of_mu(X, N)
        array_of_mu = zeros(1, N);
        for i = 1 : N
            array_of_mu(i) = get_mu(X(1 : i));
        end
    end

    function array_of_var = get_array_of_var(X, N)
        array_of_var = zeros(1, N);
        for i = 1 : N
            array_of_var(i) = get_Ssqr(X(1 : i));
        end
    end
    
    function mu_higher = get_mu_higher(array_of_mu, array_of_var, alpha, N)
        mu_higher = zeros(1, N);
        for i = 1 : N
            mu_higher(i)  = array_of_mu(i) + sqrt(array_of_var(i) ./ i) .* tinv(1 - alpha, i - 1);
        end
    end

    function mu_lower = get_mu_lower(array_of_mu, array_of_var, alpha, N)
        mu_lower = zeros(1, N);
        for i = 1 : N
            mu_lower(i) = array_of_mu(i) + sqrt(array_of_var(i) ./ i) .* tinv(alpha, i - 1);
        end
    end
    
    function sigma_higher = get_sigma_higher(array_of_var, alpha, N)
        sigma_higher = zeros(1, N);
        for i = 1:N
            sigma_higher(i) = array_of_var(i) .* (i - 1) ./ chi2inv(alpha, i - 1);
        end
    end
    
    function sigma_lower = get_sigma_lower(array_of_var, alpha, N)
        sigma_lower = zeros(1, N);
        for i = 1 : N
            sigma_lower(i) = array_of_var(i) .* (i - 1) ./ chi2inv(1 - alpha, i - 1);
        end
    end
    
    function sigma = get_sigma_sqr(X)
        tempMu = get_mu(X);
        sigma = sum((X - tempMu) .* (X - tempMu)) / size(X,2);
    end
    
end