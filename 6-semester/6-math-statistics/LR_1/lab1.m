function lab1
    clear
    
    X = [3.38,1.21,1.85,2.24,4.17,2.99,4.81,2.71,2.70,4.41,3.21,3.15,2.77,4.05,3.89,1.56,2.78,2.04,2.82,3.28,2.63,1.89,3.57,3.15,3.80,5.40,3.25,2.04,2.61,5.06,2.87,2.66,4.80,3.86,0.09,2.45,2.40,2.14,1.69,2.36,5.44,2.77,1.94,2.55,3.97,1.88,3.01,4.21,4.74,2.02,2.38,2.46,3.51,2.89,1.57,3.53,0.77,3.31,3.58,2.77,3.61,3.71,2.38,3.06,4.29,4.76,1.69,1.59,3.21,2.74,3.99,3.53,3.52,2.84,1.21,2.82,4.34,3.65,2.22,2.87,3.14,3.58,1.96,3.41,3.85,1.96,3.02,4.22,3.10,2.68,3.67,1.70,5.47,5.02,2.52,3.09,2.19,4.44,2.33,2.27,3.34,3.05,4.35,3.58,3.43,4.49,3.57,3.20,1.53,3.53,3.53,1.27,3.40,4.53,2.21,3.28,3.50,2.01,3.30,1.86];
   
    Mmin = min(X)
    Mmax = max(X)
    R = Mmax - Mmin
    mu = getmu(X)
    Ssqr = getS(X)
    m = getNIntervals(size(X, 2))

    % graph1
    findIntervals(X, m)
    hold on;
    f(X, mu, Ssqr, mu, R);
    legend('bar chart','density function');
    hold off;

    % graph2
    figure;
    empericalF(sort(X));
    hold on;
    F(sort(X), mu, Ssqr, m, R);
    grid on;
    legend('empirical distribution function','distribution function');
    legend('Location','southeast')
    hold off;

    function mu = getmu(X)
        mu = sum(X)/size(X,2);
    end

    function sigma = getSigmaSqr(X)
        tempMu = getmu(X);
        sigma = sum((X - tempMu) .* (X - tempMu))/size(X,2);
    end

    function Ssqr = getS(X)
        n = size(X,2);
        Ssqr = n / (n - 1) * getSigmaSqr(X);
    end

    function m = getNIntervals(size)
        m = floor(log2(size)) + 2;
    end

    function findIntervals(X, m)
        sortX = sort(X);
        n = size(sortX,2);
        delta = (sortX(end) - sortX(1)) / m;
        J = sortX(1):delta:sortX(end)
        nEl = zeros(1, m);

        for i = 1:n
            for j = 1:(size(J,2) - 1)
                if (sortX(i) >= J(j) && sortX(i) < J(j+1))
                    nEl(j) = nEl(j) + 1;
                    break;
                end
            end
        end
        nEl(end) = nEl(end) + 1

        for i = 1:size(nEl,2)
            nEl(i) = nEl(i)/(n * delta);
        end
        J = [J(1) J];
        nEl = [0 nEl 0];

        stairs(J, nEl), grid;
    end

    function f(X, MX, DX, m, R)
        delta = R/m;
        sigma = sqrt(DX);

        Xn = min(X):delta/20:max(X);
        Y = normpdf(Xn, MX, sigma);
        plot(Xn, Y, '-.');
    end

    function empericalF(X)  
        [yy, xx] = ecdf(X);
        stairs(xx, yy);
    end

    function F(X, MX, DX, m, R)
        delta = R/m;
        Xn = min(X):delta/20:max(X);
        Y = 1/2 * (1 + erf((Xn - MX) / sqrt(2*DX))); 

        plot(Xn, Y, '--');
    end
end