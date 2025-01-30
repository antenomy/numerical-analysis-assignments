%% Assignment a)
x = -10:.001:10;

y = x.^2 - (8 * x) - 10 * sin( (3.5 * x) + 1) + 20;

% plot(x, y)

tolerance = 1e-6;
root_indices = find(abs(y) < tolerance);

roots_x = x(root_indices);

disp('Roots of y = x^2:');
disp(roots_x);

% Hur många nollställen finns det? 6 stycken
% x_1 ≈ 1.81
% x_2 ≈ 2.17
% x_3 ≈ 3.44
% x_4 ≈ 4.08
% x_5 ≈ 5.27 
% x_6 ≈ 5.77

%% Assignment b)
format long

function y = func1(x)
    y = x.^2 - (8 * x) - 10 * sin( (3.5 * x) + 1) + 20;
end

function xn = fixpunkt(x, iteration)
    xn = 0.05 * (x.^2 + (12*x) - (10 * sin((3.5 * x) + 1 ))) + 1;
    
    if abs(func1(xn)) > 10^(-10)
        iteration = iteration + 1;
        xn = fixpunkt(xn, iteration);
    else
        disp(['Iteration: ', num2str(iteration, '%.11g'), '    xn: ', num2str(xn, '%.11g'), '    y: ', num2str(func1(xn), '%.11g')])
    end
end

num_to_check = [1.81, 2.17, 3.44, 4.08, 5.27, 5.77];

for num = num_to_check
    fixpunkt(num, 1);
end

% 1.81  Iteration: 8     xn: 1.8152602476    y:  1.8310686301 e-11
% 2.17  Iteration: 17    xn: 1.8152602476    y:  3.6084912836 e-11
% 3.44  Iteration: 48    xn: 3.4323780522    y: -8.3229423353 e-11
% 4.08  Iteration: 58    xn: 3.4323780522    y:  8.2469142626 e-11
% 5.27  Iteration: 16    xn: 5.2701564597    y:  5.3269388900 e-11
% 5.77  Iteration: 28    xn: 5.2701564597    y: -5.9436899846 e-11

%% Assignment c)



%% Assignment d)


