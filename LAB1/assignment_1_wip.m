%% Assignment a)
x = -10:.001:10;

y = x.^2 - (8 * x) - 10 * sin( (3.5 * x) + 1) + 20;

plot(x, y)

tolerance = 1e-6;
root_indices = find(abs(y) < tolerance);

roots_x = x(root_indices);

disp('Roots of y = x^2:');
disp(roots_x);

% Question 1: How many roots are there?
% Answer: There are 6. They are:
% x_1 ≈ 1.81
% x_2 ≈ 2.17
% x_3 ≈ 3.44
% x_4 ≈ 4.08
% x_5 ≈ 5.27 
% x_6 ≈ 5.77



%% Assignment b)
format long

global print_list iteration_count
print_list = zeros(1, 100);
iteration_count = 0;

function y = func1(x)
    y = x.^2 - (8 * x) - 10 * sin( (3.5 * x) + 1) + 20;
end

function x_n = fixpoint(x, original_x, iteration, last_ten)
    global print_list iteration_count
    x_n = 0.05 * (x.^2 + (12*x) - (10 * sin((3.5 * x) + 1 ))) + 1;

    if last_ten == 1
        print_list(iteration) = abs(x_n - x);
        iteration_count = iteration;
    end
       
    
    if abs(func1(x_n)) > 10^(-10)
        iteration = iteration + 1;
        x_n = fixpoint(x_n, original_x, iteration, last_ten);
    else
        disp(['Starting Point: ', num2str(original_x), '   Total Iterations: ', num2str(iteration, '%.11g'), '    Approximate root: ', num2str(x_n, '%.11g')])
        if last_ten == 1
            iteration_count = iteration;
        end
    end
end

num_to_check = [1.81, 2.17, 3.44, 4.08, 5.27];

for num = num_to_check
    fixpoint(num, num, 1, 0);
end

fixpoint(5.77, 5.77, 1, 1);

fprintf('\nLast 10 |xn+1 - xn| values for starting point 5.77:\n')

print_list = print_list(iteration_count-9:iteration_count);

for element = print_list
    disp(element)
end




% Program Output:
% Starting Point: 1.81   Total Iterations: 8    Approximate root: 1.8152602476
% Starting Point: 2.17   Total Iterations: 17    Approximate root: 1.8152602476
% Starting Point: 3.44   Total Iterations: 48    Approximate root: 3.4323780522
% Starting Point: 4.08   Total Iterations: 58    Approximate root: 3.4323780522
% Starting Point: 5.27   Total Iterations: 16    Approximate root: 5.2701564597
% Starting Point: 5.77   Total Iterations: 28    Approximate root: 5.2701564597
% 
% Last 10 |xn+1 - xn| values for starting point 5.77:
%     2.537049956785609e-07
%     8.150794172934184e-08
%     2.618612526816833e-08
%     8.412836649540623e-09
%     2.702798873599477e-09
%     8.683302965550865e-10
%     2.789688480220320e-10
%     8.962430797510024e-11
%     2.879296800983866e-11
%     9.250378241176804e-12

% Question 1: Why can this new equation be used as a fix point iteration for the previous one?
% Answer: Because the new function is the previous function rewritten as
% g(x) = x. This was done by:
%  1. Adding 20x to each side
%  2. Dividing both sides by 20

% Question 2: Which of these roots can fixed-point iteration find?
% Answer: The convergence rate for this method is equal to |g'(r)| and is
% locally convergent for |g'(r)| < 1. Thus the algorithm found the three 
% roots which fulfill this and for the three roots which didn't, the
% starting point led the algorithm to the nearest root which was
% convergent. This can be confirmed by examining the derivatives at the
% starting points:
%  x_1: ≈ -0.08
%  x_2: ≈  1.99
%  x_3: ≈ -0.64
%  x_4: ≈  2.60
%  x_5: ≈ -0.32
%  x_6: ≈  2.40
% Which correlates with the theorems assertion.

% Question 3: How does the difference |x_n+1 - x_n| decrease according to
% theory. 
% 
%
% We can verify this by 
%



%% Assignment c)

function y = func1_derivative(x)
    y = 2*x - 8  - 35 * cos( (3.5 * x) + 1);
end

function x_n = newton(x, original_x, iteration, print_diff)
    x_n = x - (func1(x)/func1_derivative(x));

    if print_diff == 1
        disp(abs(x_n - x))
    end
    
    if abs(func1(x_n)) > 10^(-10)
        iteration = iteration + 1;
        x_n = newton(x_n, original_x, iteration, print_diff);
    else
        disp(['Starting Point: ', num2str(original_x), '   Total Iterations: ', num2str(iteration, '%.11g'), '    Approximate root: ', num2str(x_n, '%.11g')])
    end
    
end

num_to_check_2 = [1.81, 2.17, 3.44, 5.27, 5.77];

for num = num_to_check_2
    newton(num, num, 1, 0);
end

fprintf('\n|xn+1 - xn| values for starting point 5.77:\n')

newton(4.08, 4.08, 1, 1);

% Program Output:
% Starting Point: 1.81   Total Iterations: 3    Approximate root: 1.8152602476
% Starting Point: 2.17   Total Iterations: 3    Approximate root: 2.1714433738
% Starting Point: 3.44   Total Iterations: 3    Approximate root: 3.4323780522
% Starting Point: 5.27   Total Iterations: 2    Approximate root: 5.2701564597
% Starting Point: 5.77   Total Iterations: 2    Approximate root: 5.7704911928
%
% |xn+1 - xn| values for starting point 4.08:
%     0.004492833875932
%     1.635993122395263e-05
%     2.120641440228610e-10
%
% Starting Point: 4.08   Total Iterations: 3    Approximate root: 4.0844764737



%% Assignment d)

format long

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

newton_array = zeros(1,80);
newton_array(1) = 1.81;
iteration_n = 0;
x_n = 1.81;
e_n = 1;

while e_n > 10^(-16)
    x_n = newton_d(x_n);
    e_n = abs(x_n - reference_x);

    newton_array(iteration_n+2) = e_n;
    iteration_n = iteration_n + 1;
end

fixpoint_array = zeros(1,80);
fixpoint_array(1) = 1.81;
iteration_f = 0;
x_n = 1.81;
e_n = 1;

while e_n > 10^(-16)
    x_n = fixpoint_d(x_n);
    e_n = abs(x_n - reference_x);

    fixpoint_array(iteration_f+2) = e_n;
    iteration_f = iteration_f + 1;
end

max_iteration = max(iteration_f, iteration_n);



n = 1:1:max_iteration+1;



e_m = newton_array(x);
e_f = fixpoint_array(x);
% Maybe we should turn this funtion into one that directly uses the amount
% of iterations through netwons method

semilogy(n, e_m)
hold on
semilogy(n, e_f)

% Reference Value of x: 1.815260247632966