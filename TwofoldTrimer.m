%{
    --------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [February 2026]
    Description:  [Tridiagonal 3-Toeplitz with complex entries]
    --------------------------------------------------------------
%}

clc
clear all;
close all;

% --- Parameters ---
    n = 30;  % size of matrix (should be multiple of three)

% --- Build alternating diagonal values ---  
    a1 = 1.2 + 0.0*1i;
    a2 = 1.0 + 0.0*1i;
    a3 = 1.4 + 0.1*1i;

    b1 = 3.8 + 0.4*1i;
    b2 = 1.2 - 0.9*1i;
    b3 = 2.0 + 0.1*1i;
    c1 = b1;
    c2 = b2; 
    c3 = b3;


    B0 = [a2, b2;
          c2, a3];
    sigma_B0 = eig(B0);

% --- Generate 2-Toeplitz matrix ---
    main_diag = repmat([a1; a2; a3], n/3, 1);
    
    % --- Above diagonal ---
    above_diag = repmat([b1; b2; b3], n/3 -1, 1);
    if length(above_diag) < n-1
        above_diag(end+1) = b1;
        above_diag(end+1) = b2;
    end
    
    % --- Below diagonal ---
    below_diag = repmat([c1; c2; c3], n/3 - 1, 1);
    if length(below_diag) < n-1
        below_diag(end+1) = c1;
        below_diag(end+1) = c2;
    end
    
    T = diag(main_diag) + diag(above_diag, 1) + diag(below_diag, -1);

    % --- Generate the submatrix for twofold structure ---
    T_Sub   = T(2:n, 2:n);
    eigvals = eig(T_Sub);

% --- Generate chiral symmetric Twofold Toeplitz matrix ---

    % --- Generate the reflection ---
    T_flip = rot90(T,2).';
    
    % --- Polulate Twofold operator ---
    T_mirror = zeros(2*n-1);
    T_mirror(1:n, 1:n) = T_flip;
    T_mirror(n:end, n:end) = T;
    
    eta = -3*T(1,1);
    q   = T(1,2);
    T_mirror(n,n) = eta; 

    % --- Compute the Spectrum ---
    eigMirror = eig(T_mirror);

    [VMir, DiaMir] = eig(T_mirror);


    % --- Generate the submatrix ---
    EigSub = eig(T_Sub);

    % --- Impedence matched Gap Frequency ---

    % Initial guess (important!)
    lambda0 = -20 + 1i;
    x0 = [real(lambda0); imag(lambda0)];

    lambda01 = 10 + 1i;
    x01 = [real(lambda01); imag(lambda01)];
    
    % Options
    opts = optimoptions('fsolve','Display','iter','TolFun',1e-12,'TolX',1e-12);
    
    % Solve
    [xsol,  ~, ~] = fsolve(@(x) fsolve_g_wrapper(x,q,eta,a1,a2,a3,b1,b2,b3,c1,c2,c3), x0,  opts);
    [xsol2, ~, ~] = fsolve(@(x) fsolve_g_wrapper(x,q,eta,a1,a2,a3,b1,b2,b3,c1,c2,c3), x01, opts);
    % Recover complex lambda
    lambda_sol   = xsol(1)  + 1i*xsol(2);
    lambda_sol_2 = xsol2(1) + 1i*xsol2(2);
    lambda_sol = [lambda_sol; lambda_sol_2];
    fprintf('Solved lambda = %.15f + %.15fi\n', real(lambda_sol), imag(lambda_sol));

    % --- Plot the spectrum ---
    figure;
    plot(nan, nan, 'bx', 'MarkerSize', 8, 'LineWidth', 4.5);
    hold on;
    plot(nan, nan, 'ro', 'MarkerSize', 8, 'LineWidth', 4.5);
    plot(real(lambda_sol), imag(lambda_sol), 'ko', 'MarkerSize', 12, 'MarkerFaceColor', 'g');
    plot(real(sigma_B0), imag(sigma_B0), 'ks', 'MarkerSize', 17, 'MarkerFaceColor', 'c');

    %plot(real(lambda_roots(1,:)), imag(lambda_roots(1,:)), 'k-', 'LineWidth', 2);
    %plot(real(lambda_roots(2,:)), imag(lambda_roots(2,:)), 'k-', 'LineWidth', 2);
    plot(real(eigMirror), imag(eigMirror), 'bx', 'MarkerSize', 8, 'LineWidth', 1.5);
    plot(real(EigSub),    imag(EigSub),    'ro', 'MarkerSize', 8, 'LineWidth', 1.5);
    xlabel('$\mathrm{Re}$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('$\mathrm{Im}$', 'Interpreter', 'latex', 'FontSize', 14);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18);
    set(gcf, 'Position', [100, 100, 500, 300]); 

    xlim([-1 + min( real(eigMirror) ), 1 + max( real(eigMirror) )]);
    ylim([-0.5 + min( imag(eigMirror) ), 0.7 + max( imag(eigMirror) )]);
    grid on;
    legend({'$\sigma(\mathbf{T}_{AB})$','$\sigma(\mathbf{T}_{A})$','$\lambda_{int}$','$\sigma(\mathbf{B}_0)$'}, ...
       'Interpreter','latex','Location','northeast','NumColumns',4);
    hold off;
 
%% --- Defining functions ---
  
function gval = g(lambda,q,eta,a1,a2,a3,b1,b2,b3,c1,c2,c3)

    % --- Roots via discriminat ---
    a  = c1*c2*c3;
    A0 = [a1-lambda, b1, 0; c1, a2-lambda, b2; 0, c2, a3-lambda];
    b  = det(A0) - b3*c3*(a2-lambda);
    c  = b1*b2*b3;
    
    D = b^2 - 4*a*c;          % Discriminant
    
    z1 = (-b + sqrt(D)) / (2*a);
    z2 = (-b - sqrt(D)) / (2*a);
    
    z_all = [z1; z2];

    % ---- Step 2: Pick root of smallest magnitude ----
    [~, idx] = min(abs(z_all));
    z_star   = z_all(idx);

    M = [a1,          b1, c3 / z_star;
         c1,          a2, b2;
         b3 * z_star, c2, a3       ];

    %disp(det(M-lambda*eye(3)));

    [~,~,V] = svd(M-lambda*eye(3));
    v = V(:,end);   % singular vector for smallest singular value
    ratio = v(1)/v(3);

    disp(v)
    eigvec = (M-lambda*eye(3))*v;
    disp(eigvec);

    gval = 2*q*z_star*ratio + eta;
  
end


function Fvec = fsolve_g_wrapper(x,q,eta,a1,a2,a3,b1,b2,b3,c1,c2,c3)
    % x(1) = Re(lambda), x(2) = Im(lambda)
    lambda = x(1) + 1i*x(2);
    
    F = g(lambda,q,eta,a2,a3,a1,b2,b3,b1,c2,c3,c1) - lambda;
    
    % fsolve expects real vector, so split real/imag
    Fvec = [real(F); imag(F)];
end
