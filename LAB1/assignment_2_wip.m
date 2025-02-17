

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

function allSolutions=findSolutions(ra, rb, a, b)
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
end

function partB()
    ra=0.8;
    rb=1.2;
    a=[-1; 2.5];
    b=[2; 0];

    
    
    allSolutions = findSolutions(ra, rb, a, b);
    
    for j = 1:size(allSolutions, 2)
        fprintf('Solution %d:\n', j);
        disp(allSolutions(:, j));
    end
    

    plotSolutions(allSolutions, a, b, ra, rb);
end
%X0 = [0;0;0;0];
%ra=0.8;
%rb=1.2;
%a=[-1; 2.5];
%b=[2; 0];

%[Xrot, iter] = punkter(X0, ra, rb, a, b);
%disp('Solution:');
%disp(Xrot);
%disp(['Iterations: ', num2str(iter)]);
%partB()

%%%%%%%%% PLOT:





%%%%%%%%%%%% Part (c)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%plot(xTriangle, yTriangle, '-o', 'LineWidth', 2);
%hold on;
%fill(xTriangle, yTriangle, 'cyan', 'FaceAlpha', 0.3);
%hold off;
%axis equal; 
%grid on;

function hasX=checkPolyX(T1,T2,C1,C2,C3,r1,r2,r3)
    % T1 and T2 are the end points of the tangent line solution
    % C1,C2,C3 are centerpoints of circle as usual
    xTriangle = [C1(1),C2(1),C3(1),C1(1)];
    yTriangle = [C1(2),C2(2),C3(2),C1(2)];
    xTangent = [T1(1),T2(1)];
    yTangent = [T1(2),T2(2)];
    [xi,yi] = polyxpoly(xTriangle,yTriangle,xTangent,yTangent);

    % Plot the triangle
    plotTriangle = false;
    if plotTriangle == true
        figure;
        hold on;
        axis equal;
        grid on;
        xlabel('X');
        ylabel('Y');
        title('intersect check');
        plot(xTangent, yTangent, 'r-', 'LineWidth', 2);
        plot(xTriangle, yTriangle, '-o', 'LineWidth', 2);
        rectangle('Position', [C1(1)-r1, C1(2)-r1, 2*r1, 2*r1], 'Curvature', [1,1], 'EdgeColor', 'w', 'LineWidth', 1);
        rectangle('Position', [C2(1)-r2, C2(2)-r2, 2*r2, 2*r2], 'Curvature', [1,1], 'EdgeColor', 'w', 'LineWidth', 1);
        rectangle('Position', [C3(1)-r3, C3(2)-r3, 2*r3, 2*r3], 'Curvature', [1,1], 'EdgeColor', 'w', 'LineWidth', 1);
      
        if ~isempty(xi)
            plot(xi, yi, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'yellow'); % Highlight intersections
        end
        hold off;
    end
    

    hasX=~isempty(xi);
end

function outerSolution=getOuterSolution(r1,r2,rOuter,p1,p2,pOuter)
    allSolutions=findSolutions(r1,r2,p1,p2);
    for i = 1:size(allSolutions, 2)
        solution = allSolutions(:, i);
        T1 = solution(1:2);
        T2 = solution(3:4);
        if ~checkPolyX(T1,T2,p1,p2,pOuter,r1,r2,rOuter);
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

    % Angles of arc midpoints between contact points wrt circle
    %thetaM1 = (thetaP1 + thetaP2) / 2;
    %thetaM2 = thetaM1 + pi;
    % Midpoint coords
    %M1 = [C(1) + rC * cos(thetaM1), C(2) + rC * sin(thetaM1)];
    %M2 = [C(1) + rC * cos(thetaM2), C(2) + rC * sin(thetaM2)];
    % secant lines between first contact point and midpoints
    %secantP1M1 = [M1(1)-P1(1), M1(2)-P1(2)];
    %secantP1M2 = [M2(1)-P1(1), M2(2)-P1(2)];
    % secant line between contact points
    %secantP1P2 = P2-P1;
    % 2D "cross" product of midpoint secants with contact secant 
    %crossM1 = secantP1P2(1)*secantP1M1(2)-secantP1P2(2)*secantP1M1(1);
    %crossM2 = secantP1P2(1)*secantP1M2(2)-secantP1P2(2)*secantP1M2(1);
    
    %if crossM1 > 0 % If 1st midpoint is on outer side
    %    arcStart = thetaP1;
    %    arcEnd = thetaP2;
    %else % Otherwise 2nd midpoint must be outer side
    %    arcStart = thetaP2;
    %    arcEnd = thetaP1;
    %end
    %arcDiff = mod(arcEnd - arcStart, 2*pi); 

    %if startTheta > endTheta
    %    arcSubArray1 = linspace(startTheta, pi, 30);
    %    arcSubArray2 = linspace(-pi, endTheta, 30);
    %    arcArray = [arcSubArray1, arcSubArray2];
    %else
    %    arcArray = linspace(startTheta, endTheta, 50);
    %end
    %xArc = C(1) + rC * cos(arcArray);
    %yArc = C(2) + rC * sin(arcArray);
    %outerArcDist = abs(startTheta - endTheta); %CCW from start_theta to end_theta
    %disp(outerArcDist);
end


function L=solveStringSystem(ABsol,BCsol,CAsol,ra,rb,rc,a,b,c)
    AB1=ABsol(1:2); AB2=ABsol(3:4);
    BC1=BCsol(1:2); BC2=BCsol(3:4);
    CA1=CAsol(1:2); CA2=CAsol(3:4);
    [arcDistA, xArcA, yArcA] = getOuterArcDist(CA2,AB1,ra,a);
    [arcDistB, xArcB, yArcB] = getOuterArcDist(AB2,BC1,rb,b);
    [arcDistC, xArcC, yArcC] = getOuterArcDist(BC2,CA1,rc,c);

    ABDist = norm([AB1(1) AB2(1)]-[AB1(2) AB2(2)]);
    BCDist = norm([BC1(1) BC2(1)]-[BC1(2) BC2(2)]);
    CADist = norm([CA1(1) CA2(1)]-[CA1(2) CA2(2)]);

    L = ABDist + BCDist + CADist + arcDistA + arcDistB + arcDistC;
    plotStringSystem = true;
    if plotStringSystem == true
        
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
        %disp(arcDistA);
        
        % Plot the two arc midpoint candidates
        
        hold off;
    end
end

function [L,X]=langd(ra,rb,rc,a,b,c)
    % Indata:
    %
    % ra - radie för cirkel a (skalär)
    % rb - radie för cirkel a (skalär)
    % rc - radie för cirkel a (skalär)
    % a - kolumnvektor med koordinater för mittpunkt a (xa,ya)^T
    % b - kolumnvektor med koordinater för mittpunkt b (xb,yb)^T
    % c - kolumnvektor med koordinater för mittpunkt c (xc,yc)^T
    %
    % Utdata:
    %
    % L - längden på snöret (skalär)
    % X - De tre lösningsvektorerna samlade i en
    % matris (4 rader, 3 kolumner), X=[Xab, Xbc, Xca]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % kod...

    ABsol = getOuterSolution(ra,rb,rc,a,b,c);
    BCsol = getOuterSolution(rb,rc,ra,b,c,a);
    CAsol = getOuterSolution(rc,ra,rb,c,a,b);

    L = solveStringSystem(ABsol,BCsol,CAsol,ra,rb,rc,a,b,c);
    X = [ABsol, BCsol, CAsol];
end

ra = 0.6;
rb = 0.9;
rc = 1.5;
a = [-1; 2];
b = [3; 0.5];
c = [0.4; -2.5];
[L,X]=langd(ra,rb,rc,a,b,c);
disp(L);
disp(X);