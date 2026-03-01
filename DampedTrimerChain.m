%{
    --------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [February 2026]
    Description:  [Open resonator chain with complex wavespeeds]
    --------------------------------------------------------------
%}

clc
clear;
close all;

% --- Number of resonators ---
    n = 60;

% --- Resonator spacings ---
    s1 = 1.0;
    s2 = 1.5;
    s3 = 1.2;

% --- Wavespeeds ---
    v1 = 1 - 1.0*1i;
    v2 = 1 - 1.0*1i;
    v3 = 1 - 1.0*1i;

% --- Generate the Capacitance matrix ---
    C = generate_Capacitance(n, s1, s2, s3);

% --- Generalised Capacitance ---
    d = repmat([v1, v2, v3], 1, n/3);
    V = diag(d); 
    G = V * C;

    EigG = eig(G);

% --- Generate the open limit ---
    
    % --- Extract the coefficients ---
        a1 = G(4,4);
        a2 = G(2,2);
        a3 = G(3,3);
        b1 = G(1,2);
        b2 = G(2,3);
        b3 = G(3,4);
        c1 = G(2,1);
        c2 = G(3,2);
        c3 = G(4,3);

        recip = sqrt( (b1 * b2 * b3) / (c1 * c2 * c3) );

        N = 100;
        alpha_vals   = linspace(0, 2*pi, N);
        lambda_roots = zeros(3, N);
        lambda_roots2 = zeros(3, N);

        % --- Solve polynomial in lambda ---
        for k = 1:N
                alphaVal = alpha_vals(k);
        
                z_star = recip * exp(-1i*alpha_vals(k));
            
                M = [a1,          b1, c3 / z_star;
                     c1,          a2, b2;
                     b3 * z_star, c2, a3];
        
                lambda_new = eig(M);
            
                if k == 1
                    lambda_roots(:,k) = lambda_new;
                else
                    lambda_old = lambda_roots(:,k-1);
            
                    % Match by minimal distance
                    cost = abs(lambda_old - lambda_new.');
                    [~, assignment] = min(cost, [], 2);
            
                    lambda_roots(:,k) = lambda_new(assignment);
                end
        
        end


% --- Plot the spectrum ---
    figure;
    plot(nan, nan, 'rx', 'MarkerSize', 8, 'LineWidth', 4.5);
    hold on;

    plot(real(lambda_roots(1,:)), imag(lambda_roots(1,:)), 'k-', 'LineWidth', 4);
    plot(real(lambda_roots(2,:)), imag(lambda_roots(2,:)), 'k-', 'LineWidth', 4);
    plot(real(lambda_roots(3,:)), imag(lambda_roots(3,:)), 'k-', 'LineWidth', 4);

    plot(real(EigG),     imag(EigG),     'rx', 'MarkerSize', 8, 'LineWidth', 1.5);

    xlabel('$\mathrm{Re}$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('$\mathrm{Im}$', 'Interpreter', 'latex', 'FontSize', 14);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18);
    set(gcf, 'Position', [100, 100, 500, 300]); 

    legend({'$\sigma(\mathcal{C}) $','$\Lambda$'}, ...
       'Interpreter','latex','Location','northeast','NumColumns',2);
    xlim([-0.2 + min( real(EigG )), 0.2 + max( real(EigG) )]);
    ylim([-0.2 + min( imag(EigG )), 0.2 + max( imag(EigG) )]);
    set(gcf, 'Position', [100, 100, 500, 300]); 
    grid on;


%% --- Defining functions ---

function C = generate_Capacitance(n, s1, s2, s3)
    % Objective: compute the interface capacitance matrix

    if mod(n, 3) ~= 0 || n <= 0
        error('Matrix size n must be of the form 3*N, where N is a natural number (N > 1).');
    end

    % Define the parameters based on the input equations
    beta1  = -1 / s1;       
    beta2  = -1 / s2;
    beta3  = -1 / s3;
    alpha1 =  1 / s1 + 1 / s2; 
    alpha2 =  1 / s2 + 1 / s3;
    alpha3 =  1 / s3 + 1 / s1;   

    C = zeros(n, n);

    % --- Populate diagonal ---
    for i = 1:n
        if mod(i, 3) == 1
            C(i, i) = alpha3;
        elseif mod(i, 3) == 2
            C(i, i) = alpha1;
        else
            C(i, i) = alpha2;
        end
    end

    % --- Populate upper and lower diagonal ---
    for i = 1:n-1
        if mod(i, 3) == 1
            C(i, i+1) = beta1; % Upper diagonal
            C(i+1, i) = beta1; % Lower diagonal       
        elseif mod(i, 3) == 2
            C(i, i+1) = beta2; % Upper diagonal
            C(i+1, i) = beta2; % Lower diagonal
        else
            C(i, i+1) = beta3; % Upper diagonal
            C(i+1, i) = beta3; % Lower diagonal
        end
    end

    % --- Manually adjust the corners ---
    
    % Manually set the corners
    C(1, 1) = - C(1, 2);
    C(n, n) = - C(n, n-1);

end

