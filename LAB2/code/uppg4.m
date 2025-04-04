%% Assignment a)


function y = func1(x)
    y = 2/(sqrt(pi))*exp(-x^2)
end

interval_start

x = 1
y = erf(x)

plot()

comparison_value = erf()



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



f = trapets(0, 1, @func1, N);



%% Assignment b)




%% Assignment c)






