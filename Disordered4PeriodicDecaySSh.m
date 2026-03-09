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
    n = 40;  % size of matrix (should be even)  

    a1 = 0;
    a2 = 0;
    a3 = 0;
    a4 = 0; 

    b1 = 2-0.4*1i;
    b2 = 0.5+0.3*1i;
    b3 = 1.3;
    b4 = 1.8-0.2*1i;

    decay = abs( (b2*b4) / (b1*b3) );

    disorder_steps = 20;
    iter = 1000;

    % --- Symmetric off diagonal entries ---
    c1 = b1; 
    c2 = b2; 
    c3 = b3;
    c4 = b4;

    % --- Common interface ---
    eta = 0;


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

% --- Generate symmetric Twofold Toeplitz matrix ---

    % --- Generate the reflection ---
    T_flip = rot90(T,2).';
    
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
    
        parfor ind = 1:iter
            T_mirror = perturb_tridiag(T_Pris, disorder);
    
            [VMir, DiaMir] = eig(T_mirror);
            [eigValsSorted, idx] = sort(diag(DiaMir), 'ascend');
            VMir = VMir(:, idx);
            eigMirror = eigValsSorted;
    
            [~, idx_min] = min(abs(eigMirror));
            v = VMir(:, idx_min);
    
            left_val(ind)  = -(log(abs(v(n-1))) - log(abs(v(1))))   / ((n-1)/4);
            right_val(ind) = -(log(abs(v(n+1))) - log(abs(v(end)))) / ((n-1)/4);
        end
    
        expec_vals_left(j)  = exp(mean(left_val));
        expec_vals_right(j) = exp(mean(right_val));
    
    end

    %% --- Plot the result ---

    %figure;
    %semilogy(1:length(v), abs(v), 'ko-', 'LineWidth', 4);
    %hold off;

    %%

    g = decay - 0.005 + 0 .* d_vals;

    figure;
    plot(d_vals, g, 'k-', 'LineWidth', 4);
    hold on;
    plot(d_vals, expec_vals_right, 'b-', 'LineWidth', 4);
    plot(d_vals, expec_vals_left, 'r-', 'LineWidth', 4);

    ylim([decay - 0.025, decay + 0.025]);
    xlabel('Disorder $d$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('Average decay rate', 'Interpreter', 'latex', 'FontSize', 14);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18);
    legend({'$|a/b|$','Left Average', 'Right Average'}, ...
       'Interpreter','latex','Location','northeast','NumColumns',3);
    set(gcf, 'Position', [100, 100, 500, 300]); 
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

