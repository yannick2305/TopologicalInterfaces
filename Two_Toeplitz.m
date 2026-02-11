
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
    n = 40;  % size of matrix (should be even)

    phase = 1;
    %V = diag( exp(1i*phase) * ones(n,1) );

    % Build alternating diagonal values
    d = repmat([exp(1i*phase), exp(-1i*phase/2)], 1, n/2);
    
    % Construct diagonal matrix
    V = diag(d);
    
    % --- Toeplitz matrix entries ---
    a1 = 2.2; 
    a2 = 1.2;          
    b1 = 2; 
    b2 = 2;
    c1 = 0.5*b1;
    c2 = 0.5*b2;

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

    T = V * T;

    a1 = T(3,3); 
    a2 = T(2,2);          
    b1 = T(1,2); 
    b2 = T(2,3);
    c1 = T(2,1);
    c2 = T(3,2);

    [V, Dia] = eig(T);
    Toep = T;

% --- Generate and plot the spectrum ---

    
    eigvals = eig(T);

    %
     beta = 0.5 * log( abs( (b1*b2) / (c1*c2) ));
    % Alpha range
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
        
            % Coefficients of quadratic: lambda^2 - (a1 + a2)*lambda + C
            A = 1;
            B = -(a1 + a2);
            D = B^2 - 4*A*C;
        
            % Quadratic formula
            lambda1 = (-B + sqrt(D)) / (2*A);
            lambda2 = (-B - sqrt(D)) / (2*A);
        
            lambda_roots(:, k) = [lambda1; lambda2];
        end
        
    % Plot in complex plane
    figure;
    plot(nan, nan, 'rx', 'MarkerSize', 8, 'LineWidth', 4.5);
    hold on;
    plot(real(lambda_roots(1,:)), imag(lambda_roots(1,:)), 'k-', 'LineWidth', 4);
    plot(real(lambda_roots(2,:)), imag(lambda_roots(2,:)), 'k-', 'LineWidth', 4);
      
    plot(real(eigvals), imag(eigvals), 'rx', 'LineWidth', 8, 'LineWidth', 4.5);
    xlabel('$\mathrm{Re}$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('$\mathrm{Im}$', 'Interpreter', 'latex', 'FontSize', 14);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18);
    set(gcf, 'Position', [100, 100, 500, 300]); 

    xlim([-2, 4]);
    grid on;
    axis equal;
    legend('Eigenvalues', 'Collapsed Symbol', 'Interpreter', 'latex', 'Location', 'northwest');
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


%% --- Defect modes ---

    figure;
    for d = 0:0.1:7.99
        T = Toep;
        Def = eye(n);
        defect_site = 2;
        Def(defect_site, defect_site) =  Def(defect_site, defect_site) + d;
        T = Def * T;
        eigvals = eig(T);
    
        plot(real(eigvals), imag(eigvals), 'r.', 'MarkerFaceColor', 'r');
        hold on;
    end
    
    xlabel('Real Part');
    ylabel('Imaginary Part');
    grid on;
    axis equal;
    hold off;
