
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
    n = 28;  % size of matrix (should be even)

    % --- Build constant phase shift ---
    phase = 0.0;
    V = diag( exp(1i*phase) * ones(n,1) );

    % --- Build alternating diagonal values ---
    %d = repmat([exp(1i*phase), exp(-1i*phase/2)], 1, n/2);
    %V = diag(d);   

    a1 = 1.2 + 0.0*1i;
    a2 = 1+ 0.0*1i;
    b1 = 3.8 + 0.4*1i;
    b2 = 1.2 - 0.9*1i;
    c1 = b1; %1.2 - 0.0*1i;
    c2 = b2; %0.5 + 0*1i;

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
    T = V * T;

    [V, Dia] = eig(T);

    % --- Generate the submatrix ---
    T_Sub   = T(2:n, 2:n);
    eigvals = eig(T_Sub);

% --- Generate and plot the spectrum ---

    % --- Generate the open spectrum ---
        beta = 0.5 * log( abs( (b1*b2) / (c1*c2) ));
        N = 1000;
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

    % --- Plot in complex plane ---
    figure;
    plot(nan, nan, 'rx', 'MarkerSize', 8, 'LineWidth', 4.5);
    hold on;
    plot(real(lambda_roots(1,:)), imag(lambda_roots(1,:)), 'k-', 'LineWidth', 4);
    plot(real(lambda_roots(2,:)), imag(lambda_roots(2,:)), 'k-', 'LineWidth', 4);
    plot(real(eigvals),     imag(eigvals),     'bx', 'MarkerSize', 8, 'LineWidth', 1.5);
    xlabel('$\mathrm{Re}$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('$\mathrm{Im}$', 'Interpreter', 'latex', 'FontSize', 14);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18);
    set(gcf, 'Position', [100, 100, 500, 300]); 

    xlim([-0.2 + min(min(real(lambda_roots(1,:))), min(real(lambda_roots(2,:)))), 0.2 + max(max(real(lambda_roots(1,:))), max(real(lambda_roots(2,:))))]);
    ylim([-3.8, 3.8]);
    set(gcf, 'Position', [100, 100, 500, 300]); 
    grid on;
    
    legend('Eigenvalues', 'Collapsed Symbol', 'Interpreter', 'latex', 'Location', 'northwest');
    hold off;

%% --- Generate chiral symmetric Twofold Toeplitz matrix ---

    % --- Generate the reflection ---
    T_flip = rot90(T,2).';
    
    % --- Polulate Twofold operator ---
    T_mirror = zeros(2*n-1);
    T_mirror(1:n, 1:n) = T_flip;
    T_mirror(n:end, n:end) = T;
    
    eta = 2*T(1,1);
    q = T(1,2);
    Tap = T(2,3);
    gapFreq = -2*q^2/0.7448 + eta;
    T_mirror(n,n) = eta;  % Interface unaffected by this entry

    % --- Compute the Spectrum ---
    eigMirror = eig(T_mirror);

    [VMir, DiaMir] = eig(T_mirror);


    % --- Generate the submatrix ---
    EigSub = eig(T_Sub);

    % --- Impedence matched Gap Frequency ---

    % Initial guess (important!)
    lambda0 = 1 +1i;
    x0 = [real(lambda0); imag(lambda0)];
    
    % Options
    opts = optimoptions('fsolve','Display','iter','TolFun',1e-12,'TolX',1e-12);
    
    % Solve
    [xsol, ~, ~] = fsolve(@(x) fsolve_g_wrapper(x,q,eta,a1,a2,b1,b2,c1,c2), x0, opts);
    
    % Recover complex lambda
    lambda_sol = xsol(1) + 1i*xsol(2);
    
    fprintf('Solved lambda = %.15f + %.15fi\n', real(lambda_sol), imag(lambda_sol));

    % --- Plot the spectrum ---
    figure;
    plot(nan, nan, 'bx', 'MarkerSize', 8, 'LineWidth', 4.5);
    hold on;
    plot(nan, nan, 'ro', 'MarkerSize', 8, 'LineWidth', 4.5);
    plot(real(lambda_sol), imag(lambda_sol), 'ko', 'MarkerSize', 12, 'MarkerFaceColor', 'g');
    plot(real(a2), imag(a2), 'ks', 'MarkerSize', 17, 'MarkerFaceColor', 'c');

    %plot(real(lambda_roots(1,:)), imag(lambda_roots(1,:)), 'k-', 'LineWidth', 2);
    %plot(real(lambda_roots(2,:)), imag(lambda_roots(2,:)), 'k-', 'LineWidth', 2);
    plot(real(eigMirror), imag(eigMirror), 'bx', 'MarkerSize', 8, 'LineWidth', 1.5);
    plot(real(EigSub),    imag(EigSub),    'ro', 'MarkerSize', 8, 'LineWidth', 1.5);
    xlabel('$\mathrm{Re}$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('$\mathrm{Im}$', 'Interpreter', 'latex', 'FontSize', 14);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18);
    set(gcf, 'Position', [100, 100, 500, 300]); 

    xlim([-1.5 + min(min(real(lambda_roots(1,:))), min(real(lambda_roots(2,:)))), 1.5 + max(max(real(lambda_roots(1,:))), max(real(lambda_roots(2,:))))]);
    ylim([-1.5, 2]);
    grid on;
    legend({'$\sigma(\mathbf{T}_{AB})$','$\sigma(\tilde{\mathbf{T}}_{A})$','$\lambda_{int}$','$\sigma(\mathbf{B}_0)$'}, ...
       'Interpreter','latex','Location','northeast','NumColumns',4);
    hold off;

%% --- Generate the generalised Brillouin zone ---
   
    GBZ = zeros(2, length(eigvals));

    r_reciproc = sqrt((b1*b2)/(c1*c2));

    theta_unit_circle    = linspace(0, 2*pi, 100);
    r_circle = r_reciproc*exp(1i*theta_unit_circle);

    for index = 1:length(eigvals)
        lambda = eigvals(index);

        delta_a = -c1*c2;
        delta_b = ((a1-lambda)*(a2-lambda) - c1*b1 - b2*c2);
        delta_c = -b1*b2;
    
        z_lambda_plus  = (-delta_b + sqrt(delta_b^2 - 4*delta_a*delta_c))/(-2*delta_a); 
        z_lambda_minus = (-delta_b - sqrt(delta_b^2 - 4*delta_a*delta_c))/(-2*delta_a);

        GBZ(:, index) = [z_lambda_plus; z_lambda_minus];

    end

    figure;
    plot(real(GBZ(1,:)), imag(GBZ(1,:)), 'ro', 'LineWidth', 2);
    hold on;
    plot(real(GBZ(2,:)), imag(GBZ(2,:)), 'ro', 'LineWidth', 2);
    plot(real(r_circle), imag(r_circle), 'k-', 'LineWidth', 2);

    xlabel('$\mathrm{Re}$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('$\mathrm{Im}$', 'Interpreter', 'latex', 'FontSize', 14);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18);
    set(gcf, 'Position', [100, 100, 500, 300]); 
    legend('$f^{-1}\big(\sigma(\mathbf{T_n})\big)$', '', '$r\mathbf{T}$', 'Interpreter', 'latex', 'Location', 'northeast');

    ylim([-1.4*abs(r_reciproc), 1.4*abs(r_reciproc)]);

    axis equal;
    grid on;
    hold off;

%% --- Numericlly compute Widom's set for edge states ---

    % --- Define Discrete search area ---
    Nx = 100; Ny = 10;
    x = linspace(0, 3, Nx);
    y = linspace(-0.5, 0.5, Ny);
    [X, Y] = meshgrid(x, y);
    lambda_grid = X + 1i*Y;

    detF_grid = zeros(Ny, Nx);

    % --- Loop over the complex plane ---
    for ix = 1:Nx
        for iy = 1:Ny
            g_lambda = (a1-lambda_grid(iy, ix))*(a2-lambda_grid(iy, ix)) - b1*c1 - b2*c2;
            z1 = (g_lambda - sqrt(g_lambda^2 - 4*c1*c2*b1*b2)) / (2*c1*c2);
            z2 = (g_lambda + sqrt(g_lambda^2 - 4*c1*c2*b1*b2)) / (2*c1*c2);
            rad = (abs(z1) + abs(z2))/2;
    
            detF_grid(iy, ix) = det_integral_f(f, lambda_grid(iy, ix), rad);
        end
    end

    % Plot absolute value
    figure;
    imagesc(x, y, abs(detF_grid));
    axis xy; 
    colorbar;
    xlabel('Re(\lambda)');
    ylabel('Im(\lambda)');
    clim([0 0.1]);
    title('|det F(\lambda)| in complex plane');
    
    % Plot contour where |detF| = 0 (numerically, use small tolerance)
    figure;
    contour(X,Y,abs(detF_grid), [1e-6 1e-3], 'LineWidth', 2);  % zero contour
    xlabel('Re(\lambda)');
    ylabel('Im(\lambda)');
    title('Contour in complex \lambda-plane where det F(\lambda) = 0');
    grid on;


%% --- Defining functions ---

function detF = det_integral_f(f, lambda, r)
    % detF = det_integral_f(f, lambda)
    % Computes det( (1/2pi i) * integral_T ( (f(z)-lambda*I)^(-1) dz/z ) )
    % INPUTS:
    %   f      : function handle f(z), returns 2x2 matrix
    %   lambda : scalar
    % OUTPUT:
    %   detF   : scalar value of the determinant
    
    N = 200;  % number of points for numerical integration
    theta = linspace(0, 2*pi, N+1);
    theta(end) = [];  % avoid duplication at 2pi
    
    integral_matrix = zeros(2,2);
    
    for k = 1:N
        z = r*exp(1i*theta(k));
        InvF = (f(z) - lambda*eye(2)) \ eye(2);
        integral_matrix = integral_matrix + InvF;
    end
    
    integral_matrix = integral_matrix / N;  % approximate integral / (2pi)
    detF = det(integral_matrix);
end    


function gval = g(lambda, q, eta, a1, a2, b1, b2, c1, c2)

    % --- Roots via discriminat ---
    a = c2*c1;
    b = -((a2 - lambda)*(a1 - lambda) - b2*c2 - c1*b1);
    c = b2*b1;
    
    D = b^2 - 4*a*c;          % Discriminant
    
    z1 = (-b + sqrt(D)) / (2*a);
    z2 = (-b - sqrt(D)) / (2*a);
    
    z_all = [z1; z2];

    % ---- Step 2: Pick root of smallest magnitude ----
    [~, idx] = min(abs(z_all));
    z_star = z_all(idx);

    % --- Compute eigenvector ratio ---
    m11   = a2 - lambda;
    m12   = b2 + c1/z_star;
    ratio = -m12/m11;   % this equals v1/v2

    gval = 2*q*z_star*ratio + eta;
  
end


function Fvec = fsolve_g_wrapper(x, q, eta, a1,a2,b1,b2,c1,c2)
    % x(1) = Re(lambda), x(2) = Im(lambda)
    lambda = x(1) + 1i*x(2);
    
    F = g(lambda, q, eta, a1,a2,b1,b2,c1,c2) - lambda;
    
    % fsolve expects real vector, so split real/imag
    Fvec = [real(F); imag(F)];
end
