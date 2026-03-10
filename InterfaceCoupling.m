%{
    --------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [September 2025]
    Description:  [Tridiagonal 2-Toeplitz with complex entries]
    --------------------------------------------------------------
%}

clc
clear all;
close all;

% --- Parameters ---
    n = 30;  % size of matrix (should be even) 

    a1 = 1.2 + 0.0*1i;
    a2 = 1.0 + 0.0*1i;
    b1 = 1.8 - 0.8*1i;
    b2 = 0.5 + 1.0*1i;
    c1 = b1; 
    c2 = b2; 

    % --- Coupling ---
    q1 = 1.4;
    q2 = q1;

    % --- Parity ---
    parity = 1;

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

    % --- Complexify the entries (for Capacitance) ---

    [V, Dia] = eig(T);

    % --- Generate the submatrix ---
    T_Sub   = T(2:n, 2:n);
    eigvals = eig(T_Sub);

% --- Generate and plot the spectrum ---

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

% --- Generate chiral symmetric Twofold Toeplitz matrix ---

    % --- Generate the reflection ---
    T_flip = rot90(T,2).';
    
    % --- Polulate Twofold operator ---
    T_mirror = zeros(2*n);
    T_mirror(1:n, 1:n) = T_flip;
    T_mirror(n+1:end, n+1:end) = T;

    % --- Set Coupling ---
    T_mirror(n+1,n) = q1;
    T_mirror(n,n+1) = q1;
    
    % --- Compute the Spectrum ---
    eigMirror = eig(T_mirror);

    [VMir, DiaMir] = eig(T_mirror);


    % --- Generate the submatrix ---
    EigSub = eig(T_Sub);

    % --- Impedence matched Gap Frequency ---

    % Initial guess (important!)
    lambda0 = 2.3 + 0.01*1i;
    x0 = [real(lambda0); imag(lambda0)];

    lambda1 = -3 - 1*1i;
    x1 = [real(lambda1); imag(lambda1)];
    % Options
    opts = optimoptions('fsolve','Display','off','TolFun',1e-12,'TolX',1e-12);
    
    % Solve
    [xsol, ~, ~]  = fsolve(@(x) fsolve_g_wrapper(x,q1, 1,a1,a2,b1,b2,c1,c2), x0, opts);
    [xsolp, ~, ~] = fsolve(@(x) fsolve_g_wrapper(x,q1,-1,a1,a2,b1,b2,c1,c2), x1, opts);
    % Recover complex lambda
    lambda_sol_plus = xsol(1) + 1i*xsol(2);
    lambda_sol_neg = xsolp(1) + 1i*xsolp(2);
    
    fprintf('Solved lambda for a = +1 : %.15f + %.15fi\n', real(lambda_sol_plus), imag(lambda_sol_plus));
    fprintf('Solved lambda for a = -1 : %.15f + %.15fi\n', real(lambda_sol_neg), imag(lambda_sol_neg));

    % --- Plot the spectrum ---
    figure;
    plot(nan, nan, 'bx', 'MarkerSize', 8, 'LineWidth', 4.5);
    hold on;

    plot(real(lambda_sol_plus), imag(lambda_sol_plus), 'ko', 'MarkerSize', 12, 'MarkerFaceColor', 'g');
    plot(real(lambda_sol_neg),  imag(lambda_sol_neg),  'kd', 'MarkerSize', 12, 'MarkerFaceColor', 'g');

    plot(real(a2), imag(a2), 'ks', 'MarkerSize', 17, 'MarkerFaceColor', 'c');

    plot(real(eigMirror), imag(eigMirror), 'bx', 'MarkerSize', 8, 'LineWidth', 1.5);
     
    %plot(real(lambda_roots(1,:)), imag(lambda_roots(1,:)), 'k-', 'LineWidth', 2);
    %plot(real(lambda_roots(2,:)), imag(lambda_roots(2,:)), 'k-', 'LineWidth', 2);

    xlabel('$\mathrm{Re}$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('$\mathrm{Im}$', 'Interpreter', 'latex', 'FontSize', 14);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18);
    set(gcf, 'Position', [100, 100, 500, 300]); 

    xlim([-1.5 + min(min(real(lambda_roots(1,:))), min(real(lambda_roots(2,:)))), 1.5 + max(max(real(lambda_roots(1,:))), max(real(lambda_roots(2,:))))]);
    ylim([-2.2, 2.8]);
    grid on;
    legend({'$\sigma(\mathbf{T}_{AB})$', '$\lambda_{int}^+$', '$\lambda_{int}^-$', '$\sigma(\mathbf{B}_0)$'}, ...
       'Interpreter','latex','Location','northeast','NumColumns',4);
    hold off;

    %%

    lambda = 0.29055;
    g(lambda, q1, parity, a1, a2, b1, b2, c1, c2)

%% --- Defining functions ---

function gval = g(lambda, q, parity, a1, a2, b1, b2, c1, c2)

    % --- Roots via discriminat ---
    a = c2*c1;
    b = -((a2 - lambda)*(a1 - lambda) - b2*c2 - c1*b1);
    c = b2*b1;
    
    D = b^2 - 4*a*c;   % Discriminant
    
    z1 = (-b + sqrt(D)) / (2*a);
    z2 = (-b - sqrt(D)) / (2*a);
    
    z_all = [z1; z2];

    % ---- Step 2: Pick root of smallest magnitude ----
    [~, idx] = min(abs(z_all));
    z_star = z_all(idx);

    % --- Compute eigenvector ratio ---
    m11   = a2 - lambda;
    m12   = b2 + c1/z_star;
    ratio = -m11/m12 / z_star;   % this equals v1/v2

    disp(ratio)
    gval = (1/ratio) * b1 + parity * q + a1;

    %gval = a1 + q/parity + b1*q*ratio;
  
end


function Fvec = fsolve_g_wrapper(x, q, a, a1,a2,b1,b2,c1,c2)
    % x(1) = Re(lambda), x(2) = Im(lambda)
    lambda = x(1) + 1i*x(2);
    
    F = g(lambda, q, a, a1,a2,b1,b2,c1,c2) - lambda;
    
    % fsolve expects real vector, so split real/imag
    Fvec = [real(F); imag(F)];
end
