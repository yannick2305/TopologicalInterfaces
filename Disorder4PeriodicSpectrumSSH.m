%{
    --------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [March 2026]
    Description:  [Aperiodic 4-periodic SSH]
    --------------------------------------------------------------
%}

clc
clear all;
close all;

% --- Parameters ---
    n = 40;  % size of matrix (should be even)  

    a1 = 0;
    a2 = 0;
    a3 = 0;
    a4 = 0; 

    b1 = 2-0.4*1i;
    b2 = 0.5+0.3*1i;
    b3 = 1.3;
    b4 = 1.8-0.2*1i;

    % --- Symmetric off diagonal entries ---
    c1 = b1; 
    c2 = b2; 
    c3 = b3;
    c4 = b4;

    % --- Common interface ---
    eta = 0;

    % --- Perturbation size ---
    PertuSize = 0.3;

    decay = abs( (b2*b4) / (b1*b3) );

% --- Generate 2-Toeplitz matrix ---
    main_diag = repmat([a1; a2; a3; a4], n/4, 1);
    
    % --- Above diagonal ---
    above_diag = repmat([b1; b2; b3; b4], n/4, 1);
    if length(above_diag) == n
        above_diag(end) = [];
    end
    
    % --- Below diagonal ---
    below_diag = repmat([c1; c2; c3; c4], n/4, 1);
    if length(below_diag) == n
        below_diag(end) = [];
    end
    
    T = diag(main_diag) + diag(above_diag, 1) + diag(below_diag, -1);

    [V, Dia] = eig(T);

% --- Generate the essential spectrum ---

    N = 100;
    alpha_vals   = linspace(0, 2*pi, N);
    lambda_roots = zeros(4, N);

    % --- Solve polynomial in lambda ---
    for k = 1:N
            alphaVal = alpha_vals(k);
        
            z_star = exp(-1i*alpha_vals(k));
            
            M = [a1,          b1,  0, c4 / z_star;
                 c1,          a2, b2, 0;
                  0,          c2, a3, b3;
                  b4 * z_star, 0, c3, a4];
        
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


% --- Generate symmetric Twofold Toeplitz matrix ---

    % --- Generate the reflection ---
    T_flip = rot90(T,2).';
    
    % --- Polulate Twofold operator ---
    T_mirror = zeros(2*n-1);
    T_mirror(1:n, 1:n) = T_flip;
    T_mirror(n:end, n:end) = T;

    T_mirror(n,n) = eta; 

    % --- Compute the Spectrum ---

    % --- Perturb the twofold operator ---
    T_mirror = perturb_tridiag(T_mirror, PertuSize);

    [VMir, DiaMir] = eig(T_mirror);
    eigMirror = diag(DiaMir);
  
    % --- Plot the spectrum ---
    figure;
    plot(nan, nan, 'bx', 'MarkerSize', 8, 'LineWidth', 4.5);
    hold on;
    plot(real(a2), imag(a2), 'ks', 'MarkerSize', 17, 'MarkerFaceColor', 'c');

    plot(real(lambda_roots(1,:)), imag(lambda_roots(1,:)), 'k-', 'LineWidth', 5);
    plot(real(lambda_roots(2,:)), imag(lambda_roots(2,:)), 'k-', 'LineWidth', 5);
    plot(real(lambda_roots(3,:)), imag(lambda_roots(3,:)), 'k-', 'LineWidth', 5);
    plot(real(lambda_roots(4,:)), imag(lambda_roots(4,:)), 'k-', 'LineWidth', 5);

    plot(real(eigMirror), imag(eigMirror), 'bx', 'MarkerSize', 8, 'LineWidth', 1.5);

    xlabel('$\mathrm{Re}$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('$\mathrm{Im}$', 'Interpreter', 'latex', 'FontSize', 14);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18);
    set(gcf, 'Position', [100, 100, 500, 300]); 

    xlim([-0.5 + min(real(eigMirror)), 0.5 + max(real(eigMirror))]);
    ylim([-0.2 + min(imag(eigMirror)), 0.3 +  max(imag(eigMirror))]);
    grid on;
    legend({'$\sigma(\mathbf{T}_{AB})$','$\sigma(\mathbf{B}_0)$','$\sigma_{ess}(\mathbf{T}_{AB})$'}, ...
       'Interpreter','latex','Location','northeast','NumColumns',3);
    hold off;


%% --- Plot the interface eigenmode ---

    % --- Select interface eigenmode ---
    [~, idx_min] = min(abs(eigMirror));
    %eigvals(idx_min)
    v = VMir(:, idx_min);

    % --- Plot the eigenvector ---
    figure;
    plot(1:length(v), v, 'k-o', 'LineWidth', 2);
    hold on;
    xlabel('Index $i$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('$\mathbf{w}_i$', 'Interpreter', 'latex', 'FontSize', 14);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18);
    set(gcf, 'Position', [100, 100, 500, 300]); 
    ylim([1.2* min(real(v)), 1.2* max(real(v))]);
    grid on;

%% --- Defining functions ---  

function T_noisy = perturb_tridiag(T, d)
    % Perturb the off-diagonal entries of a tridiagonal matrix
    % T     : input tridiagonal matrix (main diagonal can be zero)
    % d     : disorder rate
    
    T_noisy = T;  
    
    [n, ~] = size(T);

    for i = 1:n-1
            r = rand;
            T_noisy(i,i+1) = T(i,i+1) * ( 1 + d*(2*r - 1) ); % above diagonal
            T_noisy(i+1,i) = T(i+1,i) * ( 1 + d*(2*r - 1) ); % below diagonal
    end
end
