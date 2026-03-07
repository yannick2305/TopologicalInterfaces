%{
    --------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [March 2026]
    Description:  [Aperiodic Dimer SSH Expected Decay rate]
    --------------------------------------------------------------
%}

clc
clear all;
close all;

% --- Parameters ---
    n = 20;                 % size of matrix (has to be even)  
    iter = 5000;            % number of realizations for average
    disorder_steps = 50;    % Discretisation of disorder 

    a1 = 0;
    a2 = 0;
    b1 = 2 + 1.0*1i;       
    b2 = 1 - 0.3*1i; 

    % --- Symmetric off diagonal entries ---
    c1 = b1; 
    c2 = b2; 

    % --- Common interface ---
    eta = 0;

% --- Generate symmetric Twofold Toeplitz matrix ---
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

    % --- Generate the reflection ---
    T_flip =  rot90(T,2).';
    
    % --- Polulate Twofold operator ---
    T_mirror = zeros(2*n-1);
    T_mirror(1:n, 1:n) = T_flip;
    T_mirror(n:end, n:end) = T;
    T_mirror(n,n) = eta; 

    T_Pris = T_mirror;

% --- Perturb the structure and measure decay ---
    d_vals = linspace(0.0001,0.99, disorder_steps);  
    expec_vals_left  = zeros(size(d_vals));
    expec_vals_right = zeros(size(d_vals));
    
    for j = 1:length(d_vals)
   
        disorder = d_vals(j);
        left_val  = zeros(1, iter);
        right_val = zeros(1, iter);

 
        for ind = 1:iter
            T_mirror = perturb_tridiag(T_Pris, b1, b2, disorder);
        
            [VMir, DiaMir] = eig(T_mirror);
            [eigValsSorted, idx] = sort(diag(DiaMir), 'ascend'); 
            VMir = VMir(:, idx);
            DiaMir = diag(eigValsSorted);
            eigMirror = diag(DiaMir);
        
            % --- Select interface eigenmode ---
            [~, idx_min] = min(abs(eigMirror));
            eigMirror(idx_min);
            v = VMir(:, idx_min);

            left_val(ind)  = -( log(abs(v(n-1))) - log( abs(v(1))))  / ((n-2)/2);
            right_val(ind) = -( log(abs(v(n+1))) - log(abs(v(end)))) / ((n-2)/2);
        end

        expec_vals_left(j)  = exp(mean(left_val));
        expec_vals_right(j) = exp(mean(right_val));

    end

    %% --- Plot the result ---

    g = abs(b2/b1) + 0 .* d_vals;

    figure;
    plot(d_vals, g, 'k-', 'LineWidth', 4);
    hold on;
    plot(d_vals, expec_vals_right, 'b-', 'LineWidth', 4);
    plot(d_vals, expec_vals_left, 'r-', 'LineWidth', 4);

    ylim([abs(b2/b1) - 0.02, abs(b2/b1) + 0.02]);
    xlabel('Disorder $d$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('Average decay rate', 'Interpreter', 'latex', 'FontSize', 14);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18);
    legend({'$|a/b|$','Left Average', 'Right Average'}, ...
       'Interpreter','latex','Location','northeast','NumColumns',3);
    set(gcf, 'Position', [100, 100, 500, 300]); 
    grid on;


%% --- Defining functions ---  


function T_noisy = perturb_tridiag(T, a, b, d)
    % Perturb the off-diagonal entries of a tridiagonal matrix
    % T     : input tridiagonal matrix (main diagonal can be zero)
    % a, b  : original off-diagonal values
    % d     : disorder rate
    
    T_noisy = T;  
    
    [n, ~] = size(T);

    for i = 1:n-1
        if T(i,i+1) == a
            r = rand; % same perturbation for above and below term
            T_noisy(i,i+1) = a * ( 1 + d*(2*r - 1) );
            T_noisy(i+1,i) = a * ( 1 + d*(2*r - 1) );
        elseif T(i,i+1) == b
            r = rand; % same perturbation for above and below term
            T_noisy(i,i+1) = b * ( 1 + d*(2*r - 1) );
            T_noisy(i+1,i) = b * ( 1 + d*(2*r - 1) );
        end
    end
end

