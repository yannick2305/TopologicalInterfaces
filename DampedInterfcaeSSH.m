%{
    --------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [March 2026]
    Description:  [SSH complex wavespeed]
    --------------------------------------------------------------
%}

clc
clear all;
close all;

% --- Set the parameters ---
    si = 10;
    n  = (4*si + 1);
    s1 = 1;
    s2 = 2;   % Set s1>s2 for edge mode.
    fs = 18;
    
% --- Generare SSH Capacitance ---
    Cap = generate_SSH_Capacitance(n, s1, s2);
    
% --- Complex wave speed ---
    phase = -1;
    V = diag( exp(1i*phase) * ones(n,1) );
    T = V * Cap;

% --- Retrieve Toeplitz structure from Capacitance ---
    a1 = T(3,3); 
    a2 = T(2,2);          
    b1 = T(2,1); 
    b2 = T(3,2);
    c1 = T(1,2);
    c2 = T(2,3);
    eta = T(floor(n/2)+1, floor(n/2)+1);

% --- Generate open limit and plot the spectrum ---

    eigvals = eig(T);
    [VecT, DiaT] = eig(T);
    [~, idx] = max(abs(eigvals));   
    eigvals(idx) = [];                 

    % --- Compute non-reciprocity ---
    beta = 0.5 * log( abs( (b1*b2) / (c1*c2) ));
    
    % --- Compute the open spectrum ---
    N = 500;
    alpha_vals = linspace(0, 2*pi, N);
    lambda_roots = zeros(2, N); % Two roots per alpha
        
    for k = 1:N
        alphaVal = alpha_vals(k);
        
        % --- Constant term of the quadratic ---
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
            
        % --- Solve quadratic equation for lambda ---
        lambda_roots(:, k) = [lambda1; lambda2];
    end
        
    % --- Ping's asymptotic interface frequency ---
    PingFreq = exp(1i*phase) * 0.5 * ( -sqrt( 9/(s1^2) - 14/(s1*s2) + 9/(s2^2) ) + 3/s1 + 3/s2 ) ;
    fprintf('Ping lambda =   %.15f + %.15fi\n', real(PingFreq), imag(PingFreq));

    % --- Impedance matched Frequency ---
    lambda0 = 0.5 - 1*1i;  % Initial guess

    % --- Solve Numerically ---
    x0 = [real(lambda0); imag(lambda0)];
    opts = optimoptions('fsolve','Display','off','TolFun',1e-12,'TolX',1e-12);
    [xsol, ~, ~] = fsolve(@(x) fsolve_g_wrapper(x,b2,eta,a1,a1,b2,b1,c2,c1), x0, opts);
    lambda_sol = xsol(1) + 1i*xsol(2);
    fprintf('Solved lambda = %.15f + %.15fi\n', real(lambda_sol), imag(lambda_sol));

    % --- Generate Figure ---
    figure;
    plot(real(lambda_roots(1,:)), imag(lambda_roots(1,:)), 'k-', 'LineWidth', 4);
    hold on;
    plot(real(lambda_roots(2,:)), imag(lambda_roots(2,:)), 'k-', 'LineWidth', 4);
    %plot(real(a2),                imag(a2),                'ks', 'MarkerSize', 17, 'MarkerFaceColor', 'c');
    plot(real(lambda_sol),        imag(lambda_sol),        'ko', 'MarkerSize', 12, 'MarkerFaceColor', 'g');
    plot(real(eigvals),           imag(eigvals),           'rx', 'MarkerSize', 10,  'LineWidth', 2.5);
    xlabel('$\mathrm{Re}$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('$\mathrm{Im}$', 'Interpreter', 'latex', 'FontSize', 14);
    ax = gca;
    ax.TickLabelInterpreter = 'latex';
    ax.FontSize = fs; 
    legend({'$\Gamma$', '', '$\lambda_{int}$', '$\sigma\bigl(\mathcal{C}_n\bigr)$'}, ...
           'Interpreter','latex','FontSize',fs,'Location','northeast','NumColumns',4);

    %$\lambda_{int}$ $\sigma(\mathbf{B_0})$
    xlim([-0.2 + min(real(eigvals)), 0.2 + max(real(eigvals))]);
    ylim([-0.3 + min(imag(eigvals)), 0.3 + max(imag(eigvals))])
    hold off;
    grid on;
    set(gcf, 'Position', [100, 100, 500, 300]);
    exportgraphics(gcf,'ComplexSSH.pdf','ContentType','vector')


%% --- Topological robustness ---

    N = 50;
    s_int_all = linspace(0.01, 5, N);
    interface_freq = zeros(1, N); % Two roots per alpha

    for i = 1:N
        s_int = s_int_all(i);
        
        q_def   = -1/s_int;
        eta_def = 2/s_int;

        mid = floor(n/2)+1;

        C_per = Cap;
        C_per(mid , mid) = eta_def;
        C_per(mid, mid +1) = q_def;
        C_per(mid ,mid -1) = q_def;
        C_per(mid +1,mid ) = q_def;
        C_per(mid -1,mid ) = q_def;

        eig_C_per = sort(eig(C_per));

        interface_freq(i) = eig_C_per(floor(n/2)+2);
    end

    % --- Generate Figure ---
    figure;
    hold on;

    xlim([0, s_int_all(end)]);
    xl = xlim;
    patch([xl(1) xl(2) xl(2) xl(1)], [1.0 1.0 2.0 2.0], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none')

    for i = 1:N
        plot(s_int_all(i), interface_freq(i), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
    end
    yline(2.0, 'r-', 'LineWidth', 1.0)
    yline(1.0, 'r-', 'LineWidth', 1.0)

    xlabel('$s_{int}$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('Interface Frequency', 'Interpreter', 'latex', 'FontSize', 14);
    ax = gca;
    ax.TickLabelInterpreter = 'latex';
    ax.FontSize = fs; 
    ylim([0.8, 2.2]);
    hold off;
    grid on;
    box on;
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
    %disp(z_star)

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
