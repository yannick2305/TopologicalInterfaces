%{
    --------------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [April 2026]
    Description:  [Connection between discrete and continuum systems]
    --------------------------------------------------------------------
%}

clear all;
clc;
close all;


% ---  Initialise resonator chain ---
    x_min = 0.0;   x_max = 1.0;
    
    C_b = 1.0;   % Background
    
    C1  = 5.0;   a1 = 0.15;   b1_inc = 0.50;   % Inclusion 1
    C2  = C1;    a2 = 0.60;   b2_inc = 0.85;   % Inclusion 2

% --- Number of discretisation points in finite difference method ---
    N = 20;

%% --- Build finite differences ---

    h = (x_max - x_min) / N;
    
    x_nodes = zeros(N+1, 1);
    for i = 0 : N
        x_nodes(i+1) = x_min + i * h;
    end
    
    x_half = zeros(N, 1);
    for i = 0 : N-1
        x_half(i+1) = x_min + (i + 0.5) * h;
    end
    
    s_half = zeros(N, 1);
    for i = 1 : N
        t = x_half(i);
        if (t >= a1) && (t <= b1_inc)
            s_half(i) = 1.0 / C1;
        elseif (t >= a2) && (t <= b2_inc)
            s_half(i) = 1.0 / C2;
        else
            s_half(i) = 1.0 / C_b;
        end
    end
    
    
    M  = N - 1;              % interior node count
    L0 = zeros(M, M);
    
    for j = 1 : M
        i = j;               % grid index of interior node j
    
        s_left  = s_half(i);       % s_{i-1/2}
        s_right = s_half(i + 1);   % s_{i+1/2}
    
        L0(j, j) = -(s_left + s_right) / h^2;     % diagonal
      
        if j > 1
            L0(j, j-1) = s_left  / h^2;           % sub-diagonal
        end
        if j < M
            L0(j, j+1) = s_right / h^2;           % super-diagonal
        end
    end
    
    % --- Coupling parameters ---
        eta = L0(1,1)/2;         
        b1  = L0(2, 1);

%% --- Floquet of finite difference matrix ---

    floquet_spectrum(L0, M)


%% --- Match the interface frequency ---

    n_re = 200;      
    n_im = 1;            
    
    % --- Create complex grid ---
    re_vals = linspace(-3.53, -2.647, n_re);   
    im_vals = linspace(0,  0, n_im);
    
    [Re_grid, Im_grid] = meshgrid(re_vals, im_vals);
    Lambda_grid = Re_grid + 1i * Im_grid;                    

% --- Compute F(lambda) over the complex grid ---

F_grid  = NaN(n_im, n_re);  
z1_grid = NaN(n_im, n_re);    

for row = 1 : n_im
    for col = 1 : n_re

        lam = Lambda_grid(row, col);

        % ── A = L_0 - lambda * I  ─────────────────────────────────────────────
        A = L0 - lam * eye(M);

        %z_pts = [1, -1, 2];
        z_pts = [1, -1, 2];
        g = arrayfun(@(z) det_fun(z, A, M) * z, z_pts);  % g(z) = z * det_fun(z)

        % g(z) = c1*z^2 + c0*z + c_neg1  --> solve 3x3 linear system
        V = [z_pts(1)^2, z_pts(1), 1;
             z_pts(2)^2, z_pts(2), 1;
             z_pts(3)^2, z_pts(3), 1];

        coeffs = V \ g(:);   % [c1; c0; c_neg1]
        c1     = coeffs(1);
        c0     = coeffs(2);
        c_neg1 = coeffs(3);
 
        
        %fprintf('Quadratic: (%.4f)z² + (%.4f)z + (%.4f) = 0\n', c1, c0, c_neg1);
        
        discriminant = c0^2 - 4*c1*c_neg1;
        sqrt_disc    = sqrt(discriminant);  
        
        if real(c0) >= 0
            q = -0.5 * (c0 + sqrt_disc);
        else
            q = -0.5 * (c0 - sqrt_disc);
        end
        
        z1 = q / c1;       
        z2 = c_neg1 / q;      


        % Double-check by substituting back into the polynomial p(z) = c1*z^2 + c0*z + c_neg1
        p = @(z) c1*z^2 + c0*z + c_neg1;
        
        %{
        fprintf('  lambda= %.2e \n', lam);

        fprintf('  z1= %.2e \n', z1);
        fprintf('  z2= %.2e \n', z2);
        fprintf('  z1*z2= %.10e \n', z2*z1);
        
        fprintf('Polynomial residuals:\n');
        fprintf('  p(z1) = %.2e + %.2ei  (should be 0)\n', real(p(z1)), imag(p(z1)));
        fprintf('  p(z2) = %.2e + %.2ei  (should be 0)\n', real(p(z2)), imag(p(z2)));
        
        % And via the original det_fun for a full end-to-end check
        fprintf('\ndet_fun residuals:\n');
        fprintf('  |det_fun(z1)| = %.2e  (should be 0)\n', abs(det_fun(z1, A, M))*z1);
        fprintf('  |det_fun(z2)| = %.2e  (should be 0)\n', abs(det_fun(z2, A, M))*z2);
        %}

        if abs(z1) < abs(z2)
            z_sol = z1;
        else
            z_sol = z2;
        end

        L_z1       = L0;                 
        L_z1(1, M) = L0(1,2) * z_sol;           
        L_z1(M, 1) = L0(1,2) / z_sol;          

        % --- Consistency check ---
    
        [~, ~, V_svd]  = svd(L_z1 - lam * eye(M));
        v              = V_svd(:, end);      % null vector  (last column of V_svd)

        v1 = v(1);     % first entry
        vk = v(M);     % last  entry  (index M = bottom of the interior vector)


        if abs(vk) < 1e-14
            continue   % avoid division by near-zero
        end

        F_grid(row, col) = (2 * L0(1,2) * z_sol * v1) / (vk)  +  eta  -  lam;

    end
end

%% --- Trace the matched interface ---

% --- Initialise the function ---
    F_re = real(F_grid);
    F_im = imag(F_grid); % should be exactly zero
    F_abs = abs(F_grid);


    A = [-min(re_vals), F_re(1)];
    B = [-max(re_vals), max(F_re)];

% --- Plot the matched frequency ---
    figure;
    plot(-re_vals, F_re, 'r-', 'LineWidth', 1.8);
    hold on;
     
    plot(A(1), A(2), 'ko', 'MarkerSize',8, 'MarkerFaceColor','k')
    plot(B(1), B(2), 'ko', 'MarkerSize',8, 'MarkerFaceColor','k')
    text(A(1)+0.02, A(2)-0.2, '  A', 'Interpreter', 'latex', 'FontSize', 16);
    text(B(1)+0.02, B(2)-0.2, '  B', 'Interpreter', 'latex', 'FontSize', 16);

    xline(-min(re_vals), '--k')   
    xline(-max(re_vals), '--k') 

    xlim([2.6, 3.6]);
    ylim([-800, 800]);
    ylabel('$F(\omega^2)$', 'Interpreter', 'latex', 'FontSize', 12);
    xlabel('$\omega^2$', 'Interpreter', 'latex', 'FontSize', 12);

    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18);
    set(gcf, 'Position', [100, 100, 500, 300]); 
    box on;
    grid on;
    hold off;


%% --- Defining functions ---
 
function floquet_spectrum(A, M)

    B0 = A(2:end,2:end);
    eigB0 = sort(real(eig(full(B0))));

    % --- Sweep theta in [0, 2*pi] ---
    N_theta = 200;                          % resolution
    theta   = linspace(-pi, pi, N_theta);
    z_vec   = exp(1i * theta);
 
    eigs_all = zeros(M, N_theta);           % store all M eigenvalues per theta
 
    for j = 1:N_theta
        Hz = bloch_matrix(A, M, z_vec(j));
        ev = eig(Hz);
        eigs_all(:, j) = sort_eigs(ev);    % sort for smooth band curves
    end
 
    % --- Plots ---
 
    A = [0, 2.648];
    B = [0, 3.554];

    figure;
    hold on;
    plot(nan, nan, 'k-',  'LineWidth', 2);
    plot(nan, nan, 'r--', 'LineWidth', 2);

    for b = 1:M
        plot(theta, -real(eigs_all(b,:)), 'k-', 'LineWidth', 2);
    end

    yline(-eigB0, 'r--', 'LineWidth', 2);

    yline(A(2), 'k:', 'LineWidth', 2);
    yline(B(2), 'k:', 'LineWidth', 2);

    plot(A(1), A(2), 'ko', 'MarkerSize',8, 'MarkerFaceColor','k')
    plot(B(1), B(2), 'ko', 'MarkerSize',8, 'MarkerFaceColor','k')
    text(A(1), A(2)-0.6, '  B', 'Interpreter', 'latex', 'FontSize', 16);
    text(B(1), B(2)+0.5, '  A', 'Interpreter', 'latex', 'FontSize', 16);

    xlabel('$\theta$', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('$\omega^2$',      'Interpreter', 'latex', 'FontSize', 12);
    legend({'$\sigma\big(f(e^{i\theta})\big)$', '$\sigma(\mathbf{B}_0)$'}, 'Interpreter','latex','Location','northeast','NumColumns',2);
    grid on;  
    box on;

    xlim([-pi-0.2, pi+0.2]);
    
    xticks(-pi:pi/2:pi)
    xticklabels({'$-\pi$','$-\pi/2$','$0$','$\pi/2$','$\pi$'});
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18);
    set(gcf, 'Position', [100, 100, 500, 300]); 

    ylim([0, 10]);
end

function Hz = bloch_matrix(A, M, z)
    Hz = A;
    Hz(1, M) = A(1,2) * z;  
    Hz(M, 1) = A(1,2) / z;  
end
 
function ev_sorted = sort_eigs(ev)
    [~, idx] = sortrows([real(ev), imag(ev)]);
    ev_sorted = ev(idx);
end

function f = det_fun(z, A, M)
    Lz_minus_lam = A;
    Lz_minus_lam(1, M) = A(1,2) * z;
    Lz_minus_lam(M, 1) = A(1,2) / z;

    f = det(Lz_minus_lam);
end

