
ra=0.8;
rb=1.2;
a=[-1; 2.5];
b=[2; 0];

function Fx = F(x, ra, rb, a, b)
    % x(1): x_1
    % x(2): y_1
    % x(3): x_2
    % x(4): y_2
    % a(1): x_a
    % a(2): y_a
    % b(1): x_b
    % b(2): y_b
    % ra: r_a
    % rb: r_b
    Fx = [(x(1)-a(1))^2+(x(2)-a(2))^2-ra^2;
        (x(3)-b(1))^2+(x(4)-b(2))^2-rb^2;
        (x(1)-x(3))*(x(1)-a(1))+(x(2)-x(4))*(x(2)-a(2));
        (x(1)-x(3))*(x(3)-b(1))+(x(2)-x(4))*(x(4)-b(2))];
    %disp(Fx);
end

function DF=jacobian(x, a, b)
    DF=[2*(x(1)-a(1)), 2*(x(2)-a(2)), 0, 0;
        0, 0, 2*(x(3)-b(1)), 2*(x(4)-b(2));
        2*x(1)-x(3)-a(1), 2*x(2)-x(4)-a(2), -x(1)+a(1), -x(2)+a(2);
        x(3)-b(1), x(4)-b(2), -2*x(3)+x(1)+b(1), -2*x(4)+x(2)+b(2)];
end



function [Xrot, iter] = punkter(X0, ra, rb, a, b)
    % Indata:
    %
    % X0 - kolumnvektor med startgissningarna (x1,y1,x2,y2)^T
    % ra - radie för cirkel a (skalär)
    % rb - radie för cirkel a (skalär)
    % a - kolumnvektor med koordinater för mittpunkt a (xa,ya)^T
    % b - kolumnvektor med koordinater för mittpunkt b (xb,yb)^T
    %
    % Utdata:
    %
    % Xrot - kolumnvektor med lösningen (x1,y1,x2,y2)^T
    % iter - antal iterationer som använts
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    TOL = 1e-10;
    Xrot = X0;
    iter = 0;
    max_iter = 100;
    while max(abs(F(Xrot, ra, rb, a, b))) > TOL && iter < max_iter
        iter = iter + 1;
        s = jacobian(Xrot, a, b) \ -F(Xrot, ra, rb, a, b);
        Xrot = Xrot + s;
    end
end

%X0 = [0;0;0;0];
%[Xrot, iter] = punkter(X0, ra, rb, a, b);
%disp('Solution:');
%disp(Xrot);
%disp(['Iterations: ', num2str(iter)]);

monte_carlo_attempts = 100;
xMin = min(a(1), b(1)) - 2*max(ra, rb);
xMax = max(a(1), b(1)) + 2*max(ra, rb);
yMin = min(a(2), b(2)) - 2*max(ra, rb);
yMax = max(a(2), b(2)) + 2*max(ra, rb);

monte_carlo_guesses = zeros(4, monte_carlo_attempts);
for k = 1:monte_carlo_attempts
    monte_carlo_guesses(1,k) = xMin + (xMax - xMin)*rand; % x_1
    monte_carlo_guesses(2,k) = yMin + (yMax - yMin)*rand; % y_1
    monte_carlo_guesses(3,k) = xMin + (xMax - xMin)*rand; % x_2
    monte_carlo_guesses(4,k) = yMin + (yMax - yMin)*rand; % y_2
end


allSolutions = [];
allIters = [];

for k = 1:monte_carlo_attempts
    X0 = monte_carlo_guesses(:, k);
    [Xrot, iter] = punkter(X0, ra, rb, a, b);
    
    % Check if it converged
    if max(abs(F(Xrot, ra, rb, a, b))) < 1e-10
        % Remove duplicates
        isDuplicate = false;
        if ~isempty(allSolutions)
            for i = 1:size(allSolutions, 2)
                difference = allSolutions(:, i) - Xrot;

                if max(abs(difference)) < 1e-7
                    isDuplicate = true;
                    break
                end
            end
        end
        if isDuplicate == false
            allSolutions = [allSolutions, Xrot];
        end
    end
end

for j = 1:size(allSolutions, 2)
    fprintf('Solution %d:\n', j);
    disp(allSolutions(:, j));
end

function plotSolutions(solutions, a, b, ra, rb)
    figure;
    hold on;
    axis equal;
    xlim([-6, 6]);
    ylim([-6, 6]);
    grid on;
    title('Solution plot');
    xlabel('x');
    ylabel('y');
    scatter(a(1), a(2), 40, 'r', 'filled');
    scatter(b(1), b(2), 40, 'b', 'filled');
    rectangle('Position', [a(1)-ra, a(2)-ra, 2*ra, 2*ra], 'Curvature', [1,1], 'EdgeColor', 'w', 'LineWidth', 1);
    rectangle('Position', [b(1)-rb, b(2)-rb, 2*rb, 2*rb], 'Curvature', [1,1], 'EdgeColor', 'w', 'LineWidth', 1);
    
    for i = 1:size(solutions, 2)
        X1 = solutions(1:2, i);
        X2 = solutions(3:4, i);
        scatter(X1(1), X1(2), 20, 'g', 'filled');
        scatter(X2(1), X2(2), 20, 'g', 'filled');
        plot([X1(1), X2(1)], [X1(2), X2(2)], 'g', 'LineWidth', 1);
    end
    
    hold off;

    %plotX = linspace(-6, 6)
    %plotY = 
end

plotSolutions(allSolutions, a, b, ra, rb);
