%{
    --------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [September 2025]
    Description:  [Robustness Complex Dimer edge mode]
    --------------------------------------------------------------
%}

clc
clear all;
close all;

% --- Parameters ---
    n = 80;  % size of matrix (should be even)
    
% --- Toeplitz matrix entries ---
    b1 = 1.0 * exp(1*1i); 
    b2 = b1*exp(-1*1i) + 0.2;

    a1 = -b1-b2; 
    a2 = a1+0.0;          
    c1 = 1*b1;
    c2 = 1*b2;

% --- Generate 2-Toeplitz matrix ---
    main_diag = repmat([a1; a2], n/2, 1);
    
    % --- Above diagonal ---
    above_diag = repmat([b1; b2], floor((n-1)/2), 1);
    if length(above_diag) < n-1
        above_diag(end+1) = b1;
    end
    
    % --- Below diagonal ---
    below_diag = repmat([c1; c2], floor((n-1)/2), 1);
    if length(below_diag) < n-1
        below_diag(end+1) = c1;
    end
    
    T = diag(main_diag) + diag(above_diag, 1) + diag(below_diag, -1);

    [V, Dia] = eig(T);
    Toep = T;

% --- Generate and plot the spectrum ---
    eigvals = eig(T);

    beta = 0.5 * log( abs( (b1*b2) / (c1*c2) ));
    N = 1000;
    alpha_vals = linspace(0, pi, N);
    lambda_roots = zeros(2, N); % Two roots per alpha
        
    for k = 1:N
        alphaVal = alpha_vals(k);
        
        % Constant term of the quadratic
        C = a1*a2 ...
            - c1*c2 * exp(-1i * alphaVal + beta) ...
            - b1*b2 * exp( 1i * alphaVal - beta) ...
            - c1*b1 - b2*c2;
        
        % --- Polynomial coefficients ---
        A = 1;
        B = -(a1 + a2);
        D = B^2 - 4*A*C;
        
        % Quadratic formula
        lambda1 = (-B + sqrt(D)) / (2*A);
        lambda2 = (-B - sqrt(D)) / (2*A);
        
        lambda_roots(:, k) = [lambda1; lambda2];
    end
        
    % --- Plot the spectrum ---
    figure;
    %plot(nan, nan, 'rx', 'MarkerSize', 8, 'LineWidth', 4.5);
    hold on;
    plot(real(lambda_roots(1,:)), imag(lambda_roots(1,:)), 'k-', 'LineWidth', 4);
    plot(real(lambda_roots(2,:)), imag(lambda_roots(2,:)), 'k-', 'LineWidth', 4);

    plot(real(a2), imag(a2), 'ks', 'MarkerSize', 17, 'MarkerFaceColor', 'c');
      
    %plot(real(eigvals), imag(eigvals), 'rx', 'LineWidth', 8, 'LineWidth', 4.5);
    xlabel('$\mathrm{Re}$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('$\mathrm{Im}$', 'Interpreter', 'latex', 'FontSize', 14);

       legend({'$\sigma_{ess}\big(\mathbf{T}(f)\big)$','', '$\sigma(\mathbf{B}_0)$'}, 'Interpreter','latex','Location','northeast','NumColumns',2);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18);
    set(gcf, 'Position', [100, 100, 500, 300]); 

    xlim([-3.5, 0.5]);
    ylim([-2, 0.5]);
    grid on;
    box on;
    hold off;
