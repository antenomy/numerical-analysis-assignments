%% Assignment a)

figure(1)

% g(x) function
function y = g_function(x)
    y = 2/(sqrt(pi))*exp(-x.^2);
end

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


N_step = 50;
N_range = N_step:N_step:500;
x_values = [0.11, 0.32, 1.14];

error_array = zeros(3, length(N_range));


% Error calculation for each N, and each x
for i = 1:length(x_values)
    for j = 1:length(error_array)
        error_array(j, i) = abs(erf(x_values(i)) - trapets(0, x_values(i), @g_function, j*N_step));
    end
end

% Plotting Code
loglog(N_range, error_array(:, 1));
xlabel('Number of Subintervals (N)');
ylabel('Error p');
grid on;
hold on;
loglog(N_range, error_array(:, 2));
hold on;
loglog(N_range, error_array(:, 3));
hold off;
G

% x = 0.11 : C = 0.1813049167
% x = 0.32 : C = 0.1349920833
% x = 1.14 : C = 0.08199691667


%% Assignment b)




%% Assignment c)






