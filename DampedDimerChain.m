%{
    -------------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [March 2026]
    Description:  [Open dimer resonator chain with complex wavespeeds]
    -------------------------------------------------------------------
%}

clear;
close all;

% --- Number of resonators ---
    n = 200;

% --- Resonator spacings ---
    s1 = 2.0;
    s2 = 1.5;

% --- Wavespeeds ---
    v1 = 1 + 1*1i;
    v2 = 1 + 1*1i;

% --- Generate the Capacitance matrix ---
    C = generate_Capacitance(n, s1, s2);

% --- Generalised Capacitance ---
    d = repmat([v1, v2], 1, n/2);
    V = diag(d); 
    G = V * C;

    EigG = eig(G);

    % --- Generate and plot the spectrum ---

    a1 = G(3,3);
    a2 = G(2,2);
    b1 = G(1,2);
    b2 = G(2,3);
    c1 = G(2,1); 
    c2 = G(3,2); 

    % --- Generate the open spectrum ---
        beta = 0.5 * log( abs( (b1*b2) / (c1*c2) ));
        N = 100;
        alpha_vals = linspace(pi, 2*pi, N);
        lambda_roots = zeros(2, N); % Two roots per alpha
        
        for k = 1:N
            alphaVal = alpha_vals(k);
        
            % Constant term of the quadratic
            C = a1*a2 ...
                - c1*c2 * exp(-1i * alphaVal + beta) ...
                - b1*b2 * exp( 1i * alphaVal - beta) ...
                - c1*b1 - b2*c2;
        
            % --- Coefficients of quadratic ---
            A = 1;
            B = -(a1 + a2);
            D = B^2 - 4*A*C;
        
            % --- Quadratic formula ---
            lambda1 = (-B + sqrt(D)) / (2*A);
            lambda2 = (-B - sqrt(D)) / (2*A);
        
            lambda_roots(:, k) = [lambda1; lambda2];
        end

% --- Plot the spectrum ---
    figure;
    plot(nan, nan, 'ro', 'MarkerSize', 8, 'LineWidth', 4.5);
    hold on;
    plot(real(lambda_roots(1,:)), imag(lambda_roots(1,:)), 'k-', 'LineWidth', 4);
    plot(real(lambda_roots(2,:)), imag(lambda_roots(2,:)), 'k-', 'LineWidth', 4);

    plot(real(EigG),     imag(EigG),     'ro', 'MarkerSize', 8, 'LineWidth', 1.5);

    xlabel('$\mathrm{Re}$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('$\mathrm{Im}$', 'Interpreter', 'latex', 'FontSize', 14);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18);
    set(gcf, 'Position', [100, 100, 500, 300]); 

    legend({'$\sigma(\mathcal{C}_N)$','$\Lambda$'}, ...
       'Interpreter','latex','Location','northwest','NumColumns',2);
    xlim([-0.2 + min( real(EigG )), 0.2 + max( real(EigG) )]);
    ylim([-0.2 + min( imag(EigG )), 0.2 + max( imag(EigG) )]);
    set(gcf, 'Position', [100, 100, 500, 300]); 
    grid on;


%% --- Defining functions ---

function C = generate_Capacitance(n, s1, s2)
    % Objective: compute the interface capacitance matrix

    if mod(n, 2) ~= 0 || n <= 0
        error('Matrix size n must be of the form 2*N, where N is a natural number (N > 1).');
    end

    % Define the parameters based on the input equations
    beta1  = -1 / s1;       
    beta2  = -1 / s2;
    alpha1 =  1 / s1 + 1 / s2; 
    alpha2 =  1 / s2 + 1 / s1;

    C = zeros(n, n);

    % --- Populate diagonal ---
    for i = 1:n
        if mod(i, 2) == 1
            C(i, i) = alpha2;
        else 
            C(i, i) = alpha1;
        end
    end

    % --- Populate upper and lower diagonal ---
    for i = 1:n-1
        if mod(i, 2) == 1
            C(i, i+1) = beta1; % Upper diagonal
            C(i+1, i) = beta1; % Lower diagonal       
        else
            C(i, i+1) = beta2; % Upper diagonal
            C(i+1, i) = beta2; % Lower diagonal
        end
    end

    % --- Manually adjust the corners ---
    C(1, 1) = - C(1, 2);
    C(n, n) = - C(n, n-1);

end

