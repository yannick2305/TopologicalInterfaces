%{
    --------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [March 2026]
    Description:  [Topological trimer compelx valued entries]
    --------------------------------------------------------------
%}

clc
clear all;
close all;

% --- Parameters ---
    n = 30;  % size of matrix (should be multiple of three)

% --- Generate locally symmetrix trimer symbol ---  

    a1 = 0.0 + 0.0*1i;
    a2 = 0.0 + 0.0*1i;
    a3 = 0.0 + 0.0*1i;

    b1 = 1.0 + 1.0*1i;
    b2 = 1.0 + 1.0*1i;
    b3 = b2*exp(pi/2*1i)+0.2;

    c1 = b1;
    c2 = b2; 
    c3 = b3;

% --- Trailing principal submatrix ---
    B0 = [a2, b2;
          c2, a3];

    [VB0, DiaB0] = eig(B0);
    sigma_B0 = eig(B0);

% --- Edge mode Bloch eigenvector ---
    lamB0 = sigma_B0(2);

    % --- Quadratic coefficients:  p*z^2 + q*z + r = 0 ---
    p = b1 * b2 * b3;
    q = (a1 - lamB0)*((a2 - lamB0)*(a3 - lamB0) - b2*c2) - b1*c1*(a3 - lamB0) - c3*(a2 - lamB0)*b3;
    r = c3 * c1 * c2;
    discriminant = q^2 - 4*p*r;

    z_star_1 = (-q + sqrt(discriminant)) / (2*p);
    z_star_2 = (-q - sqrt(discriminant)) / (2*p);

    F = [a1 - lamB0,    b1,         c3 / z_star_1;
         c1,            a2 - lamB0, b2;
         b3 * z_star_1, c2,         a3 - lamB0];

   [VF, DiaF] = eig(F);

% --- Generate 2-Toeplitz matrix ---
    main_diag = repmat([a1; a2; a3], n/3, 1);
    
    % --- Above diagonal ---
    above_diag = repmat([b1; b2; b3], n/3 -1, 1);
    if length(above_diag) < n-1
        above_diag(end+1) = b1;
        above_diag(end+1) = b2;
    end
    
    % --- Below diagonal ---
    below_diag = repmat([c1; c2; c3], n/3 - 1, 1);
    if length(below_diag) < n-1
        below_diag(end+1) = c1;
        below_diag(end+1) = c2;
    end
    
    T = diag(main_diag) + diag(above_diag, 1) + diag(below_diag, -1);

    EigT = eig(T);

    % --- Generate essential spectrum ---
    N = 100;
    alpha_vals   = linspace(0, pi, N);
    lambda_roots = zeros(3, N);
    lambda_roots2 = zeros(3, N);

    % --- Solve polynomial in lambda ---
    for k = 1:N
        alphaVal = alpha_vals(k);
        z_star = 1 * exp(-1i*alpha_vals(k));
            
        M = [a1,          b1, c3 / z_star;
             c1,          a2, b2;
             b3 * z_star, c2, a3];
        
        lambda_new = eig(M);
            
        if k == 1
            lambda_roots(:,k) = lambda_new;
        else
            lambda_old = lambda_roots(:,k-1);
            cost = abs(lambda_old - lambda_new.');
            [~, assignment] = min(cost, [], 2);
            lambda_roots(:,k) = lambda_new(assignment);
        end
    end

    % --- Plot the spectrum ---
    figure;
    hold on;
    plot(real(sigma_B0), imag(sigma_B0), 'ks', 'MarkerSize', 17, 'MarkerFaceColor', 'c');

    plot(real(lambda_roots(1,:)), imag(lambda_roots(1,:)), 'k-', 'LineWidth', 4);
    plot(real(lambda_roots(2,:)), imag(lambda_roots(2,:)), 'k-', 'LineWidth', 4);
    plot(real(lambda_roots(3,:)), imag(lambda_roots(3,:)), 'k-', 'LineWidth', 4);

    plot(real(sigma_B0), imag(sigma_B0), 'ks', 'MarkerSize', 17, 'MarkerFaceColor', 'c');

    plot(real(EigT), imag(EigT), 'rx', 'MarkerSize', 5, 'LineWidth', 4);
 
    xlabel('$\mathrm{Re}$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('$\mathrm{Im}$', 'Interpreter', 'latex', 'FontSize', 14);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18);
    set(gcf, 'Position', [100, 100, 500, 300]); 

    xlim([-0.5 + min( real(EigT) ), 0.5 + max( real(EigT) )]);
    ylim([-0.5 + min( imag(EigT) ), 0.7 + max( imag(EigT) )]);
    grid on;
    box on;
    legend({'$\sigma(\mathbf{B}_0)$','$\sigma_{ess}\big(\mathbf{T}(f)\big)$'}, ...
       'Interpreter','latex','Location','northwest','NumColumns',2);
    hold off;


