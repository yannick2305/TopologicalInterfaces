%{
    --------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [September 2025]
    Description:  [Tridiagonal Toeplitz with complex entries]
    --------------------------------------------------------------
%}

clear all;
close all;

% --- Parameters ---
    Np = 100;           % Plotting points symbol function
    n = 20;             % Size of the matrix
    a = 2 ;       % Main diagonal value
    b = 1 + 3*1i;       % Above diagonal
    c = 5 + 6*1i;      % Below diagonal
    fs = 18;            % Font size for plot annotation

% --- Compute the phase ---
    phase = 0.5*( angle(c) + angle(b) );

% --- Compute the modulus ---
    modulus = - 0.5 * ( log(abs(c)) + log(abs(b)) );

% --- Trace the collapsed symbol ---    
    theta = linspace(0, 2*pi, Np);
    collapsed_symbol = a - 2 * exp(-modulus) * exp(1i*phase) * cos(theta);
    %collapsed_symbol = a - 2 * b * cos(theta);


% --- Compute the spectrum numerically ---
    T = diag(a*ones(n,1)) + diag(b*ones(n-1,1),1) + diag(c*ones(n-1,1),-1);
    B = exp(-1i*phase) * T;
    [V, D] = eig(T);
    eigenvalues = diag(D);

% --- Plot the eigenvalues of a discrete system ---
    figure;
    plot(real(collapsed_symbol), imag(collapsed_symbol) , 'r-', 'LineWidth', 2);
    hold on;
    plot(real(eigenvalues),      imag(eigenvalues),       'bx', 'LineWidth', 2);
    xlabel('Real Part',     'Interpreter','latex','FontSize',fs);
    ylabel('Imaginary Part','Interpreter','latex','FontSize',fs);
    ax = gca;
    ax.TickLabelInterpreter = 'latex';
    ax.FontSize = fs; 
    hold off;
    grid on;
    xlim([min(real(eigenvalues))-1, max(real(eigenvalues))+1])
    ylim([min(imag(eigenvalues))-1, max(imag(eigenvalues))+1])
    set(gcf, 'Position', [100, 100, 500, 300]);



%% --- Defect Modes ---


    figure;
    for d = 0:0.1:20
        T = B;
        Def = eye(n);
        defect_site = 2;
        Def(defect_site, defect_site) =  Def(defect_site, defect_site) + d;
        %T(20,20) = T(20,20) + d;
        T = Def * T;
        eigvals = eig(T);
    
        plot(real(eigvals), imag(eigvals), 'r.', 'MarkerFaceColor', 'r');
        hold on;
    end

% --- Plot the collapsed symbol ---
    t = linspace(-10, 10, 100);  
    line_points = t * exp(1i*phase) + a;    
    plot(real(collapsed_symbol), imag(collapsed_symbol) , 'b-', 'LineWidth', 2);

    
    xlabel('Real Part');
    ylabel('Imaginary Part');
    grid on;
    axis equal;
    hold off;
