%{
    --------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [September 2025]
    Description:  [SSH complex wavespeed]
    --------------------------------------------------------------
%}

clear all;
close all;

    % --- Set the parameters ---
    si = 5;
    n = (4*si + 1);
    s1 = 2;
    s2 = 1;  % Set s1>s2 for edge mode.
    fs = 18;
    
    % --- Generare SSH Capacitance ---
    C = generate_SSH_Capacitance(n, s1, s2);
    
    % --- Complex wave speed ---
    phase = -1;
    V = diag( exp(1i*phase) * ones(n,1) );
    T = V * C;

    % --- Retrieve Toeplitz structure from Capacitance ---
    a1 = T(3,3); 
    a2 = T(2,2);          
    b1 = T(2,1); 
    b2 = T(3,2);
    c1 = T(1,2);
    c2 = T(2,3);

    % --- Retrieve Edge mode ---
    Edge_Toeplitz = T(2,2);

% --- Generate and plot the spectrum ---

    eigvals = eig(T);
    [~, idx] = max(abs(eigvals));   
    eigvals(idx) = [];                 

    %
    % --- Compute non-reciprocity ---
    beta = 0.5 * log( abs( (b1*b2) / (c1*c2) ));
    % --- Compute the open spectrum ---
        N = 500;
        alpha_vals = linspace(0, 2*pi, N);
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
            
            % --- Solve quadratic equation for lambda ---
            lambda_roots(:, k) = [lambda1; lambda2];
        end
        
    % --- Generate Figure ---
    figure;
    plot(real(lambda_roots(1,:)), imag(lambda_roots(1,:)), 'r-', 'LineWidth', 2);
    hold on;
    plot(real(lambda_roots(2,:)), imag(lambda_roots(2,:)), 'r-', 'LineWidth', 2);
    plot(real(Edge_Toeplitz), imag(Edge_Toeplitz), 'ko', 'MarkerSize', 12, 'MarkerFaceColor', 'g');
    plot(real(eigvals), imag(eigvals), 'bx', 'LineWidth', 2);
    
    xlabel('Real Part',     'Interpreter','latex','FontSize',fs);
    ylabel('Imaginary Part','Interpreter','latex','FontSize',fs);
    ax = gca;
    ax.TickLabelInterpreter = 'latex';
    ax.FontSize = fs; 
    legend({'$\bigcup_{\alpha \in [0, 2\pi)} \bigl(\tilde{f}(e^{-\mathrm{i}\alpha})\bigr)$', '', '$\sigma(\mathbf{B_0})$', '$\sigma\bigl(\mathcal{C}_n\bigr)$'}, ...
           'Interpreter','latex','FontSize',fs, 'Location', 'best');
    hold off;
    grid on;
    set(gcf, 'Position', [100, 100, 500, 300]);


%% --- Defining function ---


function C = generate_SSH_Capacitance(n, s1, s2)
    % Objective: compute the interface capacitance matrix

    if mod(n - 1, 4) ~= 0 || n <= 0
        error('Matrix size n must be of the form 4*N + 1, where N is a natural number (N > 1).');
    end

    % Define the parameters based on the input equations
    beta1 = -1 / s1;          
    beta2 = -1 / s2;         
    alpha = 1 / s1 + 1 / s2;  
    eta = 2 / s2;             
    tilde_alpha = 1 / s1;     

    C = zeros(n, n);

    % Populate diagonal
    for i = 1:n
        C(i, i) = alpha;
    end

    % Populate upper and lower diagonal before the defect
    for i = 1:floor(n/2)
        if mod(i, 2) == 1
            C(i, i+1) = beta1; % Upper diagonal
            C(i+1, i) = beta1; % Lower diagonal
        else
            C(i, i+1) = beta2; % Upper diagonal
            C(i+1, i) = beta2; % Lower diagonal
        end
    end
    
    % Populate upper and lower diagonal after the defect
    for i = floor(n/2)+1:n-1
        if mod(i, 2) == 1
            C(i, i+1) = beta2; % Upper diagonal
            C(i+1, i) = beta2; % Lower diagonal
        else
            C(i, i+1) = beta1; % Upper diagonal
            C(i+1, i) = beta1; % Lower diagonal
        end
    end   

    % Manually set the corners
    C(1, 1) = tilde_alpha;
    C(n, n) = tilde_alpha;

    % Set the defect in the middle
    if mod(n, 2) == 1
        mid = ceil(n / 2);
        C(mid, mid) = eta;
    else
        mid1 = n / 2;
        mid2 = mid1 + 1;
        C(mid1, mid1) = eta;
        C(mid2, mid2) = eta;
    end
end
