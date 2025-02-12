%% 1a - plotta f(x)

figure(1)

x = -10:.001:10;
y = x.^2 - (8 * x) - 10 * sin( (3.5 * x) + 1) + 20;

plot(x, y);

%% 1b - fixpunktiterationer

format long

global print_list iteration_count
print_list = zeros(1, 100);
iteration_count = 0;

function y = func1(x)
    y = x.^2 - (8 * x) - 10 * sin( (3.5 * x) + 1) + 20;
end

function xn = fixpoint(x, original_x, iteration, last_ten)
    global print_list iteration_count
    xn = 0.05 * (x.^2 + (12*x) - (10 * sin((3.5 * x) + 1 ))) + 1;

    if last_ten == 1
        print_list(iteration) = abs(xn - x);
        iteration_count = iteration;
    end
       
    
    if abs(func1(xn)) > 10^(-10)
        iteration = iteration + 1;
        xn = fixpoint(xn, original_x, iteration, last_ten);
    else
        disp(['Starting Point: ', num2str(original_x), '   Total Iterations: ', num2str(iteration, '%.11g'), '    Approximate root: ', num2str(xn, '%.11g')])
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

%% 1c - Newton

function y = func1_derivative(x)
    y = 2*x - 8  - 35 * cos( (3.5 * x) + 1);
end

function xn = newton(x, original_x, iteration, print_diff)
    xn = x - (func1(x)/func1_derivative(x));

    if print_diff == 1
        disp(abs(xn - x))
    end
    
    if abs(func1(xn)) > 10^(-10)
        iteration = iteration + 1;
        xn = newton(xn, original_x, iteration, print_diff);
    else
        disp(['Starting Point: ', num2str(original_x), '   Total Iterations: ', num2str(iteration, '%.11g'), '    Approximate root: ', num2str(xn, '%.11g')])
    end
    
end

num_to_check_2 = [1.81, 2.17, 3.44, 5.27, 5.77];

for num = num_to_check_2
    newton(num, num, 1, 0);
end

fprintf('\n|xn+1 - xn| values for starting point 4.08:\n')

newton(4.08, 4.08, 1, 1);

%% 1d - konvergensplottar

figure(2)

% kod...



function xit = fixpunkt(x0,tau)
    xn = func1(x0)

%  Indata:
%
%  x0  - startgissning (skalär)
%  tau - feltolerans (skalär)
%
%  Utdata:
%
%  xit - vektor av alla approximationer xit = [x0,x1,x2,x3,...]

if 

% kod...

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function xit = newton_d(x0,tau)

    if x_n < tau
        xit = newton();
    else
    end

%  Indata:
%
%  x0  - startgissning (skalär)
%  tau - feltolerans (skalär)
%
%  Utdata:
%
%  xit - vektor av alla approximationer xit = [x0,x1,x2,x3,...]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% kod...

end


e_n_array = newton_d(1.81, 0.5*10^(-15));