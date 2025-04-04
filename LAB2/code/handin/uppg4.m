%% Assignment a)

figure(1)

N_step = 50;
N_range = N_step:N_step:500;
x_values = [0.11, 0.32, 1.14];

% C values calculated by hand
% x = 0.11 : C = 0.1813049167
% x = 0.32 : C = 0.1349920833
% x = 1.14 : C = 0.08199691667

C_values = [0.1813049167, 0.1349920833, 0.08199691667];

error_array = zeros(3, length(N_range));
error_bounds_array = zeros(3, length(N_range));

% Error calculation for each N, and each x
for i = 1:length(x_values)
    for j = 1:length(error_array(i, :))
        error_array(i, j) = abs(erf(x_values(i)) - trapets(0, x_values(i), @g_function, N_range(j)));
        error_bounds_array(i, j) = C_values(i)*(x_values(i)^3 / N_range(j)^2);
    end
end

% Plotting Code
loglog(N_range, error_array(1, :));
xlabel('Number of Subintervals (N)');
ylabel('Error');
grid on;
hold on;
loglog(N_range, error_array(2, :));
hold on;
loglog(N_range, error_array(3, :));
hold on;
loglog(N_range, error_bounds_array(1, :), '--')
hold on;
loglog(N_range, error_bounds_array(2, :), '--')
hold on;
loglog(N_range, error_bounds_array(3, :), '--')
hold off;




%% Assignment b)

x_b = 0.5;
N_b = [50, 100, 200, 400, 800];

display_array = zeros(2, length(N_b));

for i = 1:length(N_b)
    display_array(1, i) = precision_division(0, x_b, N_b(i));
    display_array(2, i) = log2(display_array(1, i));
end

row_names = {'0.01', '0.005', '0.0025', '0.00125', '0.000625'};
table=array2table(display_array','VariableNames',{'Quotient', 'Degree of Precision'},'RowNames',row_names);
disp(table);


%% Assignment c)

figure(2);

x_step = 0.1;
x_range = 0:x_step:6;
N_values = [50, 120, 400];

error_array_2 = zeros(3, length(x_range));

% Error calculation for each x, and each N
for i = 1:length(N_values)
    for j = 1:length(error_array_2(i, :))
        error_array_2(i, j) = abs(erf(x_range(j)) - trapets(0, x_range(j), @g_function, N_values(i)));
    end
end

% Plotting Code
semilogy(x_range, error_array_2(1, :));
xlabel('x value');
ylabel('Error');
grid on;
hold on;
semilogy(x_range, error_array_2(2, :));
hold on;
semilogy(x_range, error_array_2(3, :));
hold off;



% g(x) Function
function y = g_function(x)
    y = 2/(sqrt(pi))*exp(-x.^2);
end

% Trapezoid Rule Function
function I = trapets(a, b, g, N)
    % a: Start of the interval
    % b: End of the interval
    % g: Function handle to integrate
    % N: Number of subintervals

    step_length = (b - a) / N;
    x = linspace(a, b, N + 1);
    y = g(x);

    I = step_length * (sum(y) - (y(1) + y(end)) / 2);

end

% Degree of Precision Function
function degree = precision_division(a, b, N)
    degree = (trapets(a, b, @g_function,  N) - trapets(a, b, @g_function,  2*N))/(trapets(a, b, @g_function,  2*N) - trapets(a, b, @g_function,  4*N));
end