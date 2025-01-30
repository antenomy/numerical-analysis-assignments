%% Assignment a)
x = -10:.001:10;

y = x.^2 - (8 * x) - 10 * sin( (3.5 * x) + 1) + 20;

% plot(x, y)

tolerance = 1e-6;
root_indices = find(abs(y) < tolerance);

roots_x = x(root_indices);

disp('Roots of y = x^2:');
disp(roots_x);


%% Assignment b)



%% Assignment c)



%% Assignment d)


