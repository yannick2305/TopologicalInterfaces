
%{
    --------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [March 2026]
    Description:  [Tridiagonal 2-Toeplitz with complex entries]
    --------------------------------------------------------------
%}

clc
clear all;
close all;

% --- Parameters ---
    n = 28;  % Size of matrix (should be even)

    a1 = 1.2 + 0.0*1i;
    a2 = 2.0 + 0.0*1i;
    b1 = 0.8 + 0.4*1i;
    b2 = 1.2 - 0.2*1i;
    c1 = b1; 
    c2 = b2; 

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
    T_Sub   = T(2:n, 2:n);
    eigvals = eig(T_Sub);

% --- Generate and plot the spectrum ---

    % --- Generate the essential spectrum ---
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

% --- Generate reflection symmetric Twofold Toeplitz matrix ---

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
    lambda0 = 5 - 0.8*1i;
    x0 = [real(lambda0); imag(lambda0)];

    lambda1 = 2 + 0.1*1i;
    x1 = [real(lambda1); imag(lambda1)];
    
    % Options
    opts = optimoptions('fsolve','Display','off','TolFun',1e-12,'TolX',1e-12);
    
    % Solve
    [xsol,  ~, ~] = fsolve(@(x) fsolve_g_wrapper(x,q,eta,a1,a2,b1,b2,c1,c2), x0, opts);
    [xsol2, ~, ~] = fsolve(@(x) fsolve_g_wrapper(x,q,eta,a1,a2,b1,b2,c1,c2), x1, opts);
    
    % Recover complex lambda
    lambda_sol  = xsol(1)  + 1i*xsol(2);
    lambda_sol2 = xsol2(1) + 1i*xsol2(2);
    
    fprintf('Solved lambda = %.15f + %.15fi\n', real(lambda_sol),  imag(lambda_sol));
    fprintf('Solved lambda = %.15f + %.15fi\n', real(lambda_sol2), imag(lambda_sol2));

    % --- Plot the spectrum ---
    figure;
    plot(nan, nan, 'bx', 'MarkerSize', 8, 'LineWidth', 4.5);
    hold on;
    plot(nan, nan, 'ro', 'MarkerSize', 8, 'LineWidth', 4.5);
    plot(real(lambda_sol), imag(lambda_sol), 'ko', 'MarkerSize', 12, 'MarkerFaceColor', 'g');
    plot(real(lambda_sol2), imag(lambda_sol2), 'ko', 'MarkerSize', 12, 'MarkerFaceColor', 'g');
    plot(real(a2), imag(a2), 'ks', 'MarkerSize', 17, 'MarkerFaceColor', 'c');

    %plot(real(lambda_roots(1,:)), imag(lambda_roots(1,:)), 'k-', 'LineWidth', 2);
    %plot(real(lambda_roots(2,:)), imag(lambda_roots(2,:)), 'k-', 'LineWidth', 2);
    plot(real(eigMirror), imag(eigMirror), 'bx', 'MarkerSize', 8, 'LineWidth', 1.5);
    plot(real(EigSub),    imag(EigSub),    'ro', 'MarkerSize', 8, 'LineWidth', 1.5);
    xlabel('$\mathrm{Re}$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('$\mathrm{Im}$', 'Interpreter', 'latex', 'FontSize', 14);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18);
    set(gcf, 'Position', [100, 100, 500, 300]); 

    xlim([-0.5 + min(real(EigSub)), 0.5 + max(real(EigSub))]);
    ylim([-0.5 + min(imag(EigSub)), 0.5 + max(imag(EigSub))]);
    grid on;
    legend({'$\sigma(\mathbf{T}_{AB})$','$\sigma(\mathbf{T}_{A})$','$\lambda_{int}$', '', '$\sigma(\mathbf{B}_0)$'}, ...
       'Interpreter','latex','Location','northeast','NumColumns',4);
    hold off;


%% --- Trace the impedence ---

    Re = linspace(-1.5, 5, 1000);
    Im = linspace(-1.5, 0.7, 1000);
    
    Fvals = zeros(length(Im), length(Re)); 
    
    % --- nested loop ---
    for i = 1:length(Im)
        for j = 1:length(Re)
            lambda = Re(j) + 1i*Im(i);
            Fvals(i,j) = ( g(lambda,q,eta,a1,a2,b1,b2,c1,c2) - lambda);
        end
    end

% --- Plot the contour ---
    figure;
    shading interp;

    hold on;
    plot(real(lambda_roots(1,:)), imag(lambda_roots(1,:)), 'b-', 'LineWidth', 6);
    plot(real(lambda_roots(2,:)), imag(lambda_roots(2,:)), 'b-', 'LineWidth', 6);
 
    plot(real(lambda_sol), imag(lambda_sol), 'ko', 'MarkerSize', 12, 'MarkerFaceColor', 'g');
    plot(real(lambda_sol2), imag(lambda_sol2), 'ko', 'MarkerSize', 12, 'MarkerFaceColor', 'g');
    plot(real(a2), imag(a2), 'ks', 'MarkerSize', 17, 'MarkerFaceColor', 'c');

    contour3(Re, Im, imag(Fvals), [0 0], 'k', 'LineWidth', 2)
    contour3(Re, Im, real(Fvals), [0 0], 'r', 'LineWidth', 2)

    plot(real(lambda_roots(1,:)), imag(lambda_roots(1,:)), 'b-', 'LineWidth', 6);
    plot(real(lambda_roots(2,:)), imag(lambda_roots(2,:)), 'b-', 'LineWidth', 6);

    axis xy
    xlabel('$\mathrm{Re}$', 'Interpreter', 'latex', 'FontSize', 16);
    ylabel('$\mathrm{Im}$', 'Interpreter', 'latex', 'FontSize', 16);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18);
    set(gcf, 'Position', [100, 100, 500, 300]); 

    xlim([min(Re), max(Re)]);
    ylim([min(Im), max(Im)]);
    grid on;
    box on;
    legend({'$\sigma_{ess}(\mathbf{T}_{AB})$', '','$\lambda_{int}$', '','$\sigma(\mathbf{B}_0)$', '$Im(F(\lambda))=0$','$Re(F(\lambda))=0$'}, ...
       'Interpreter','latex','Location','southwest','NumColumns',2);
    hold off;
    exportgraphics(gcf,'ContourImpedanceHermitian.pdf','ContentType','vector')


    %%

% --- complex grid ---
Re = linspace(0.5, 2.5, 200);
Im = linspace(-0.7,0.7, 200);

[ReGrid, ImGrid] = meshgrid(Re, Im);

Fvals = zeros(size(ReGrid));  % preallocate

% --- compute Fvals ---
for i = 1:length(Im)
    for j = 1:length(Re)
        lambda = Re(j) + 1i*Im(i);
        Fvals(i,j) = min(1000, g(lambda,q,eta,a1,a2,b1,b2,c1,c2) - lambda);
    end
end

% =====================================================
% Extract contour where Im(Fvals) = 0
% =====================================================

C = contourc(Re, Im, imag(Fvals), [0 0]);

% Parse contour matrix
idx = 1;
Xc = [];
Yc = [];

while idx < size(C,2)
    npts = C(2,idx);
    Xc = [Xc, C(1, idx+1:idx+npts)];
    Yc = [Yc, C(2, idx+1:idx+npts)];
    idx = idx + npts + 1;
end

% =====================================================
% Interpolate Re(Fvals) along that curve
% =====================================================

ReF_on_curve = interp2(ReGrid, ImGrid, real(Fvals), Xc, Yc);

% =====================================================
% Plot
% =====================================================

figure
hold on

% 3D curve: (Re(lambda), Im(lambda), Re(F(lambda)))
plot3(Xc, Yc, ReF_on_curve, 'r', 'LineWidth', 2)

% Optional: also show the zero contour in the plane
contour3(Re, Im, imag(Fvals), [0 0], 'k', 'LineWidth', 1.5)

z_level = 0;
plot3(real(lambda_roots(1,:)), ...
      imag(lambda_roots(1,:)), ...
      z_level * ones(size(lambda_roots(1,:))), ...
      'b-', 'LineWidth', 6);

plot3(real(lambda_roots(2,:)), ...
      imag(lambda_roots(2,:)), ...
      z_level * ones(size(lambda_roots(2,:))), ...
      'b-', 'LineWidth', 6);

grid on
xlabel('Re(\lambda)')
ylabel('Im(\lambda)')
zlabel('Re(F(\lambda))')

set(gca,'FontSize',14)
view(45,30)


% =====================================================
% Interpolate Re(F) on the curve
% =====================================================

ReF_on_curve = interp2(ReGrid, ImGrid, real(Fvals), Xc, Yc);

% =====================================================
% Create 2D projection
% =====================================================

%%
% Parameterize curve by arc length (smooth ordering)
s = [0 cumsum(sqrt(diff(Xc).^2 + diff(Yc).^2))];

    figure
    plot(s, ReF_on_curve, 'r', 'LineWidth', 2)
    xlabel('Parameter along Im(F)=0 curve')
    ylabel('Re(F)')
    grid on
    
    set(gcf, 'Position', [100, 100, 600, 400]); 

    xlabel('Curve along which $\mathrm{Im}\big(F(\lambda)\big)=0$', 'Interpreter', 'latex', 'FontSize', 16);
    ylabel('$\mathrm{Re}\big(F(\lambda)\big)$', 'Interpreter', 'latex', 'FontSize', 16);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18);
    exportgraphics(gcf,'CurveImZero.pdf','ContentType','vector')



%% --- Defining functions ---

function gval = g(lambda, q, eta, a1, a2, b1, b2, c1, c2)

    % --- Roots via discriminat ---
    a = c2*c1;
    b = -((a2 - lambda)*(a1 - lambda) - b2*c2 - c1*b1);
    c = b2*b1;
    
    D = b^2 - 4*a*c;  % Discriminant
    
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
