%% Assignment a)
x = -10:.001:10;

y = x.^2 - (8 * x) - 10 * sin( (3.5 * x) + 1) + 20;

plot(x, y)

tolerance = 1e-6;
root_indices = find(abs(y) < tolerance);

roots_x = x(root_indices);

disp('Roots of y = x^2:');
disp(roots_x);


%% Assignment b)
format long

global print_list iteration_count
print_list = zeros(1, 100);
iteration_count = 0;

function y = func1(x)
    y = x.^2 - (8 * x) - 10 * sin( (3.5 * x) + 1) + 20;
end

function x_n = fixpoint_b(x, original_x, iteration, last_ten)
    global print_list iteration_count
    x_n = 0.05 * (x.^2 + (12*x) - (10 * sin((3.5 * x) + 1 ))) + 1;

    if last_ten == 1
        print_list(iteration) = abs(x_n - x);
        iteration_count = iteration;
    end
       
    
    if abs(func1(x_n)) > 10^(-10)
        iteration = iteration + 1;
        x_n = fixpoint_b(x_n, original_x, iteration, last_ten);
    else
        disp(['Starting Point: ', num2str(original_x), '   Total Iterations: ', num2str(iteration, '%.11g'), '    Approximate root: ', num2str(x_n, '%.11g')])
        if last_ten == 1
            iteration_count = iteration;
        end
    end
end

num_to_check = [1.81, 2.17, 3.44, 4.08, 5.27];

for num = num_to_check
    fixpoint_b(num, num, 1, 0);
end

fixpoint_b(5.77, 5.77, 1, 1);

fprintf('\nLast 10 |xn+1 - xn| values for starting point 5.77:\n')

print_list = print_list(iteration_count-9:iteration_count);

for element = print_list
    disp(element)
end



%% Assignment c)

function y = func1_derivative(x)
    y = 2*x - 8  - 35 * cos( (3.5 * x) + 1);
end

function x_n = newton_c(x, original_x, iteration, print_diff)
    x_n = x - (func1(x)/func1_derivative(x));

    if print_diff == 1
        disp(abs(x_n - x))
    end
    
    if abs(func1(x_n)) > 10^(-10)
        iteration = iteration + 1;
        x_n = newton_c(x_n, original_x, iteration, print_diff);
    else
        disp(['Starting Point: ', num2str(original_x), '   Total Iterations: ', num2str(iteration, '%.11g'), '    Approximate root: ', num2str(x_n, '%.11g')])
    end
    
end

num_to_check_2 = [1.81, 2.17, 3.44, 5.27, 5.77];

for num = num_to_check_2
    newton_c(num, num, 1, 0);
end

fprintf('\n|xn+1 - xn| values for starting point 5.77:\n')

newton_c(4.08, 4.08, 1, 1);



%% Assignment d)

figure(2)

function x_n = newton_reference(x)
    x_n = x - (func1(x)/func1_derivative(x));
    
    if abs(func1(x_n)) > 10^(-16)
        x_n = newton_reference(x_n);
    end
end

function x_n = fixpoint_d(x)
    x_n = 0.05 * (x.^2 + (12*x) - (10 * sin((3.5 * x) + 1 ))) + 1;   
end

function x_n = newton_d(x)
    x_n = x - (func1(x)/func1_derivative(x));
end

reference_x = newton_reference(1.81);
% Reference Value of x: 1.815260247632966


function xit = fixpoint(x0,tau)

    max_iter = 100;
    zero_xit = zeros(1,max_iter);
    zero_xit(1) = x0;

    xn = x0;
    iter = 1;

    while abs(func1(xn)) > tau && max_iter > iter
        xn = fixpoint_d(xn);
        iter = iter + 1;
        zero_xit(iter) = xn;
    end

    xit = zero_xit(1:iter-1);

end


function xit = newton(x0,tau)

    max_iter = 100;
    zero_xit = zeros(1,max_iter);
    zero_xit(1) = x0;

    xn = x0;
    iter = 1;

    while abs(func1(xn)) > tau && max_iter > iter
        xn = newton_d(xn);
        iter = iter + 1;
        zero_xit(iter) = xn;
    end

    xit = zero_xit(1:iter-1);
end


e_newton_array = newton(1.81, 0.5*10^(-15));
e_fixpoint_array = fixpoint(1.81, 0.5*10^(-15));

semilogy(0:length(e_newton_array)-1, abs(e_newton_array-reference_x), '-o');
hold on;
semilogy(0:length(e_fixpoint_array)-1, abs(e_fixpoint_array-reference_x), '-o');

xlabel('Iteration Number');
ylabel('Error');

figure(3)

e_newton_array_n = e_newton_array(1:end-1);
e_fixpoint_array_n = e_fixpoint_array(1:end-1);

e_newton_array_n1 = e_newton_array(2:end);
e_fixpoint_array_n1 = e_fixpoint_array(2:end);

loglog(abs(e_newton_array_n-reference_x), abs(e_newton_array_n1-reference_x), '-o');
hold on;
loglog(abs(e_fixpoint_array_n-reference_x), abs(e_fixpoint_array_n1-reference_x), '-o');

xlabel('e_n');
ylabel('e_{n+1}');