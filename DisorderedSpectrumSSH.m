%{
    --------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [March 2026]
    Description:  [Aperiodic Dimer SSH]
    --------------------------------------------------------------
%}

clc
clear all;
close all;

% --- Parameters ---
    n = 30;  % size of matrix (should be even)  

    s1 = 1;
    s2 = 4;

    a1 = 0;
    a2 = 0;
    b1 = -1/s1 + 0.0*1i;
    b2 = -1/s2;

    % --- Symmetric off diagonal entries ---
    c1 = 1*b1; 
    c2 = 1*b2; 

    % --- Common interface ---
    eta = 0;

    % --- Perturbation size ---
    PertuSize = 0.6;

    % --- Check for edge states in finite system ---
    if abs(c1/b2) - abs(c2/b1) < 0
        disp('Submatrix has edge State')
    else
        disp('Submatrix has no edge State')
    end

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

    % --- Generate the submatrix ---
    T_Sub    = T(2:n, 2:n);
    eigvals  = eig(T_Sub);
    eigvalsT = eig(T);

% --- Generate and plot the spectrum ---

    % --- Generate the essential spectrum ---
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

% --- Generate symmetric Twofold Toeplitz matrix ---

    % --- Generate the reflection ---
    T_flip = rot90(T,2).';
    
    % --- Polulate Twofold operator ---
    T_mirror = zeros(2*n-1);
    T_mirror(1:n, 1:n) = T_flip;
    T_mirror(n:end, n:end) = T;
   
    q = T(1,2);
    Tap = T(2,3);
    T_mirror(n,n) = eta; 

    % --- Compute the Spectrum ---

    % --- Perturb the twofold operator ---
    Pertu = random_no_diag(2*n-1, PertuSize);
    T_mirror = T_mirror + Pertu;

    [VMir, DiaMir] = eig(T_mirror);
    eigMirror = diag(DiaMir);
  
    % --- Plot the spectrum ---
    figure;
    plot(nan, nan, 'bx', 'MarkerSize', 8, 'LineWidth', 4.5);
    hold on;
    plot(real(a2), imag(a2), 'ks', 'MarkerSize', 17, 'MarkerFaceColor', 'c');

    plot(real(lambda_roots(1,:)), imag(lambda_roots(1,:)), 'k-', 'LineWidth', 4);
    plot(real(lambda_roots(2,:)), imag(lambda_roots(2,:)), 'k-', 'LineWidth', 4);

    plot(real(eigMirror), imag(eigMirror), 'bx', 'MarkerSize', 8, 'LineWidth', 1.5);

    xlabel('$\mathrm{Re}$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('$\mathrm{Im}$', 'Interpreter', 'latex', 'FontSize', 14);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18);
    set(gcf, 'Position', [100, 100, 500, 300]); 

    xlim([-0.5 + min(real(eigMirror)), 0.5 + max(real(eigMirror))]);
    ylim([-0.2 + min(imag(eigvals)),    0.3+  max(imag(eigvals))]);
    grid on;
    legend({'$\sigma(\mathbf{T}_{AB})$','$\sigma(\mathbf{B}_0)$','$\sigma_{ess}(\mathbf{T}_{AB})$'}, ...
       'Interpreter','latex','Location','northeast','NumColumns',3);
    hold off;


%% --- Plot the interface eigenmode ---

    % --- Select interface eigenmode ---
    [~, idx_min] = min(abs(eigvals));
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

function A = random_no_diag(n, eps)
    % Generate random entries for the three diagonals
    main_diag = zeros(n,1);
    off_diag  = -eps + 2*eps*rand(n-1,1);

    % Construct tridiagonal matrix
    A = diag(main_diag) ...
      + diag(off_diag, 1) ...
      + diag(off_diag, -1);
end

