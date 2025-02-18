%% 2b)
figure(1);

ra=0.8;
rb=1.2;
a=[-1; 2.5];
b=[2; 0];

[allSolutions, allIterations] = findSolutions(ra, rb, a, b);

disp('Solutions')

for j = 1:size(allSolutions, 2)
    fprintf('Solution %d:\n', j);
    disp(allSolutions(:, j));
    disp(['Iterations: ', num2str(allIterations(j))]); 
end

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

for i = 1:size(allSolutions, 2)
    X1 = allSolutions(1:2, i);
    X2 = allSolutions(3:4, i);
    scatter(X1(1), X1(2), 20, 'g', 'filled');
    scatter(X2(1), X2(2), 20, 'g', 'filled');
    plot([X1(1), X2(1)], [X1(2), X2(2)], 'g', 'LineWidth', 1);
end

hold off;





%% 2c)

ra = 0.6;
rb = 0.9;
rc = 1.5;
a = [-1; 2];
b = [3; 0.5];
c = [0.4; -2.5];

[L,X]=langd(ra,rb,rc,a,b,c);
fprintf('Length: %.5g\n', L);
fprintf('The solution matrix:\n');
disp(X);




%% 2d)

[LDiffs, E_y, minDiff, maxDiff, iMin, iMax] = ...
    perturbationAnalysis(ra, rb, rc, a, b, c, 0.01);
%[LDiffs, forwardError, minDiff, maxDiff, iMin, iMax] = ...
%    perturbationAnalysis(ra,rb,rc,a,b,c,-0.01);

disp('Length differences:');
disp(LDiffs);

disp('E_y:');
disp(E_y);

fprintf('Smallest change = %.5g (caused by argument %d)\n', ...
    minDiff, iMin);
fprintf('Largest change = %.5g (caused by argument %d)\n', ...
    maxDiff, iMax);





%%%%% FUNCTIONS

function Fx = F(x, ra, rb, a, b)
    Fx = [(x(1)-a(1))^2+(x(2)-a(2))^2-ra^2;
        (x(3)-b(1))^2+(x(4)-b(2))^2-rb^2;
        (x(1)-x(3))*(x(1)-a(1))+(x(2)-x(4))*(x(2)-a(2));
        (x(1)-x(3))*(x(3)-b(1))+(x(2)-x(4))*(x(4)-b(2))];
end

function DF=jacobian(x, a, b)
    DF=[2*(x(1)-a(1)), 2*(x(2)-a(2)), 0, 0;
        0, 0, 2*(x(3)-b(1)), 2*(x(4)-b(2));
        2*x(1)-x(3)-a(1), 2*x(2)-x(4)-a(2), -x(1)+a(1), -x(2)+a(2);
        x(3)-b(1), x(4)-b(2), -2*x(3)+x(1)+b(1), -2*x(4)+x(2)+b(2)];
end

function [Xrot, iter] = punkter(X0, ra, rb, a, b)
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

function [allSolutions, allIterations]=findSolutions(ra, rb, a, b)
    monte_carlo_attempts = 60;
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
    allIterations = [];
    
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
                allIterations = [allIterations, iter];
            end
        end
    end
end

%%%% PART C

function hasX = checkLineX(x1,y1,x2,y2,x3,y3,x4,y4)
    hasX = false;
    denom = (y4-y3)*(x2-x1) - (x4-x3)*(y2-y1);
    if abs(denom) < eps
        return;
    end
    t = ((x4-x3)*(y1-y3) - (y4-y3)*(x1-x3)) / denom;
    u = ((x2-x1)*(y1-y3) - (y2-y1)*(x1-x3)) / denom;

    % Intersection if 0 <= t <= 1 and 0 <= u <= 1
    if (t >= 0 && t <= 1 && u >= 0 && u <= 1)
        hasX = true;
    end
end

function hasX = checkPolyX(T1,T2,C1,C2,C3)
    % T1 and T2 are the end points of the tangent line solution
    % C1,C2,C3 are centerpoints of circle as usual
    xTriangle = [C1(1),C2(1),C3(1),C1(1)];
    yTriangle = [C1(2),C2(2),C3(2),C1(2)];
    xTangent = [T1(1),T2(1)];
    yTangent = [T1(2),T2(2)];
    % Check for intersections between centerpoint triangle and tangent line
    for i = 1:3
        hasX = checkLineX( ...
            xTriangle(i),   yTriangle(i),   ...
            xTriangle(i+1), yTriangle(i+1), ...
            xTangent(1),    yTangent(1),    ...
            xTangent(2),    yTangent(2));
        if hasX 
            break;
        end
    end
end

function outerSolution=getOuterSolution(r1,r2,rOuter,p1,p2,pOuter)
    allSolutions=findSolutions(r1,r2,p1,p2);
    for i = 1:size(allSolutions, 2)
        solution = allSolutions(:, i);
        T1 = solution(1:2);
        T2 = solution(3:4);
        if ~checkPolyX(T1,T2,p1,p2,pOuter);
            outerSolution = solution;
        end
    end
end

function wrapped = wrapAngle(target, reference)
    % Fix angle range
    wrapped = target;
    while wrapped < reference
        wrapped = wrapped + 2*pi;
    end
end

function [arcDist, xArc, yArc] = getOuterArcDist(P1,P2,rC,C)
    % Angles of contact points wrt circle
    thetaP1 = atan2(P1(2)-C(2),P1(1)-C(1)); 
    thetaP2 = atan2(P2(2)-C(2),P2(1)-C(1));

    secantP1P2 = P2-P1;
    secantP1C = C-P1;
    crossValue = secantP1P2(1)*secantP1C(2) - secantP1P2(2)*secantP1C(1);

    if crossValue > 0 % If midpoint is on outer arc
        arcStart = thetaP1;
        arcEnd = thetaP2;
    else % Otherwise 2nd arc must be outer arc
        arcStart = thetaP2;
        arcEnd = thetaP1;
    end
    arcAngles = linspace(arcStart, wrapAngle(arcEnd,arcStart), 100);
    xArc = C(1) + rC * cos(arcAngles);
    yArc = C(2) + rC * sin(arcAngles);
    angleDiff = mod(arcEnd - arcStart, 2*pi);
    arcDist = rC * angleDiff;
end

function L = solveStringSystem(ABsol,BCsol,CAsol,ra,rb,rc,a,b,c)
    AB1=ABsol(1:2); AB2=ABsol(3:4);
    BC1=BCsol(1:2); BC2=BCsol(3:4);
    CA1=CAsol(1:2); CA2=CAsol(3:4);
    [arcDistA, xArcA, yArcA] = getOuterArcDist(CA2,AB1,ra,a);
    [arcDistB, xArcB, yArcB] = getOuterArcDist(AB2,BC1,rb,b);
    [arcDistC, xArcC, yArcC] = getOuterArcDist(BC2,CA1,rc,c);
       
    ABDist = norm(AB2-AB1);
    BCDist = norm(BC2-BC1);
    CADist = norm(CA2-CA1);

    L = ABDist + BCDist + CADist + arcDistA + arcDistB + arcDistC;

    plotStringSystem(ABsol,BCsol,CAsol,ra,rb,rc,a,b,c,xArcA,yArcA,xArcB,yArcB,xArcC,yArcC);
end

function plotStringSystem(ABsol,BCsol,CAsol,ra,rb,rc,a,b,c,xArcA,yArcA,xArcB,yArcB,xArcC,yArcC)
    figure;
    hold on;
    axis equal;
    grid on;
    xlabel('X');
    ylabel('Y');
    title('String System');
    ABtangentX = [ABsol(1), ABsol(3)];
    ABtangentY = [ABsol(2), ABsol(4)];
    BCtangentX = [BCsol(1), BCsol(3)];
    BCtangentY = [BCsol(2), BCsol(4)];
    CAtangentX = [CAsol(1), CAsol(3)];
    CAtangentY = [CAsol(2), CAsol(4)];
    plot(ABtangentX, ABtangentY, 'r-', 'LineWidth', 2);
    plot(BCtangentX, BCtangentY, 'r-', 'LineWidth', 2);
    plot(CAtangentX, CAtangentY, 'r-', 'LineWidth', 2);
    rectangle('Position', [a(1)-ra, a(2)-ra, 2*ra, 2*ra], 'Curvature', [1,1], 'EdgeColor', 'w', 'LineWidth', 1);
    rectangle('Position', [b(1)-rb, b(2)-rb, 2*rb, 2*rb], 'Curvature', [1,1], 'EdgeColor', 'w', 'LineWidth', 1);
    rectangle('Position', [c(1)-rc, c(2)-rc, 2*rc, 2*rc], 'Curvature', [1,1], 'EdgeColor', 'w', 'LineWidth', 1);
    text(a(1), a(2), 'a', 'FontSize', 12, 'Color', 'w');
    text(b(1), b(2), 'b', 'FontSize', 12, 'Color', 'w');
    text(c(1), c(2), 'c', 'FontSize', 12, 'Color', 'w');
    text(ABsol(1), ABsol(2), ' AB1', 'FontSize', 12, 'Color', 'g');
    text(ABsol(3), ABsol(4), ' AB2', 'FontSize', 12, 'Color', 'g');
    text(BCsol(1), BCsol(2), ' BC1', 'FontSize', 12, 'Color', 'g');
    text(BCsol(3), BCsol(4), ' BC2', 'FontSize', 12, 'Color', 'g');
    text(CAsol(1), CAsol(2), ' CA1', 'FontSize', 12, 'Color', 'g');
    text(CAsol(3), CAsol(4), ' CA2', 'FontSize', 12, 'Color', 'g');
    plot(xArcA, yArcA, 'r-', 'LineWidth', 2);
    plot(xArcB, yArcB, 'r-', 'LineWidth', 2);
    plot(xArcC, yArcC, 'r-', 'LineWidth', 2);
    hold off;
end

function [L,X]=langd(ra,rb,rc,a,b,c)

    ABsol = getOuterSolution(ra,rb,rc,a,b,c);
    BCsol = getOuterSolution(rb,rc,ra,b,c,a);
    CAsol = getOuterSolution(rc,ra,rb,c,a,b);

    L = solveStringSystem(ABsol,BCsol,CAsol,ra,rb,rc,a,b,c);
    X = [ABsol, BCsol, CAsol];
end

%%%% PART D

function [LDiffs, E_y, minDiff, maxDiff, minIndex, maxIndex] = ...
    perturbationAnalysis(ra,rb,rc,a,b,c,delta)

    lengths = zeros(1, 10);
    args = [ra, rb, rc, a(1), a(2), b(1), b(2), c(1), c(2)];
    [lengths(1), ~] = langd(ra,rb,rc,a,b,c);
    for i=1:9
        originalVal = args(i);
        args(i) = args(i) + delta;
        [lengths(i+1) ~] = langd( ...
            args(1), ...
            args(2), ...
            args(3), ...
            [args(4);args(5)], ...
            [args(6);args(7)], ...
            [args(8);args(9)] ...
        );
        args(i) = originalVal;
    end
    LDiffs = lengths(2:10) - lengths(1);
    E_y = sum(abs(LDiffs));
    [minDiff, minIndex] = min(abs(LDiffs));
    [maxDiff, maxIndex] = max(abs(LDiffs));
end