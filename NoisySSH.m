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
    
    N_dis = 100;
    disorder = linspace(0, 0.2, N_dis);
    figure;
    hold on;
    for i = 1:N_dis
        d = disorder(i);

        Cap = generate_SSH_Capacitance(n, s1, s2, d);
        eigenCap = sort(eig(Cap));

        plot(d, real(eigenCap), 'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'k');
        if s1 < s2
            plot(d, real(eigenCap(floor(n/2)+1)), 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r');
        else
            plot(d, real(eigenCap(floor(n/2)+2)), 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r');
        end

    end

    xlabel('Perturbation size $d$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('$\sigma(\mathcal{C}_n)$', 'Interpreter', 'latex', 'FontSize', 14);
    ax = gca;
    ax.TickLabelInterpreter = 'latex';
    ax.FontSize = fs; 
    ylim([-0.3, 3.5]);
    box on;
    set(gcf, 'Position', [100, 100, 500, 300]);
    hold off;

    exportgraphics(gcf,'NoisySSH.pdf','ContentType','vector');

%% --- Defining function ---

function C = generate_SSH_Capacitance(n, s1, s2, d)

  
    C = zeros(n, n);

    % Populate diagonal


    % Populate upper and lower diagonal before the defect
    for i = 1:floor(n/2)
        s2_per = s2 - d + 2*d*rand();
        s1_per = s1 - d + 2*d*rand();
       
        if mod(i, 2) == 1
            C(i, i+1) = -1/s1_per; % Upper diagonal
            C(i+1, i) = -1/s1_per; % Lower diagonal
        else
            C(i, i+1) = -1/s2_per; % Upper diagonal
            C(i+1, i) = -1/s2_per; % Lower diagonal
        end
    end
    
    % Populate upper and lower diagonal after the defect
    for i = floor(n/2)+1:n-1
        s1_per = s1 - d + 2*d*rand();
        s2_per = s2 - d + 2*d*rand();
        if mod(i, 2) == 1
            C(i, i+1) = -1/s2_per; % Upper diagonal
            C(i+1, i) = -1/s2_per; % Lower diagonal
        else
            C(i, i+1) = -1/s1_per; % Upper diagonal
            C(i+1, i) = -1/s1_per; % Lower diagonal
        end
    end   

    val = zeros(1, n-2);

    for i = 2:n-1
        C(i, i) = -C(i, i-1)-C(i, i+1);
        val(i-1) = C(i, i);
    end

    % Manually set the corners
    C(1, 1) = - C(1,2); %tilde_alpha;
    C(n, n) = -C(n-1, n); %tilde_alpha;

    % Set the defect in the middle
    mid = ceil(n / 2);
    C(mid, mid) =  -C(mid, mid-1) - C(mid, mid+1);

end
