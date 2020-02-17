%% Exercise 1: playing with correlation

wing_len = [10.4, 10.8, 11.1, 10.2, 10.3, 10.2, 10.7, 10.5, 10.8, 11.2, 10.6, 11.4];
tail_len = [7.4, 7.6, 7.9, 7.2, 7.4, 7.1, 7.4, 7.2, 7.8, 7.7, 7.8, 8.3];

figure;
plot(wing_len, tail_len, 'b*');
xlabel('Wing length (cm)');
ylabel('Tail length (cm)');

%% 1. Yes, they look positively correlated.
%% 2. Calculating Pearson's r (r_X,Y = r_Y,X since they're symmetric):

x_center = wing_len - mean(wing_len);
y_center = tail_len - mean(tail_len);

x_ssq = sum(x_center .^ 2);
y_ssq = sum(y_center .^ 2);

my_r = sum(x_center .* y_center) / sqrt(x_ssq * y_ssq);

matlab_rs = corrcoef(wing_len, tail_len);
assert(matlab_rs(1, 2) == matlab_rs(2, 1));
matlab_r = matlab_rs(1, 2);

fprintf('Using the equations, my r is %.2f; corrcoef returned %.2f.\n', my_r, matlab_r);
if my_r == matlab_r
    disp('They match!')
else
    disp('They don''t match.');
end

% They match.
%% 3. Standard error and CI

n = length(wing_len);
se = sqrt((1 - my_r^2) / (n - 2));

fprintf('The standard error of r is %.2f.\n', se);

z = 0.5 * log((1 + my_r) / (1 - my_r));
sz = sqrt(1 / (n-3));
alpha = 0.05;
z_ci = z + norminv([alpha/2, 1-(alpha/2)], 0, 1) * sz;
r_ci = (exp(2*z_ci) - 1) ./ (exp(2*z_ci) + 1);

fprintf('The 95%% confidence interval is [%.2f, %.2f].\n', r_ci(1), r_ci(2));

%% 4. Hypothesis testing

t = my_r / se;

% 2-tailed t-test: see whether r is significantly greater or less than 0.
p = 2 * (1 - tcdf(t, n-2));

disp(['The p-value is ', num2str(p), '.']);
if p < alpha
    disp('This is considered significant.');
else
    disp('This is not considered significant.');
end

%% 5. t-test with H_0: r = 0.75

r_h = 0.75;

z_m = z;
z_h = 0.5 * log((1 + r_h) / (1 - r_h));
Z = (z_m - z_h) / sqrt(1/(n - 3));
% use z-test to get p-value (?)
p = 2 * (1 - normcdf(Z));

disp(['The p-value for comparison with r=0.75 is ', num2str(p), '.']);
if p < alpha
    disp('This is considered significant.');
else
    disp('This is not considered significant.');
end

%% 6. Statistical power and sample size

r1 = 0.5;

% normalize by sem (first assuming current sample size)
se_r1 = sqrt((1 - r1^2) / (n - 2));
t1 = r1 / se_r1;

% p=0.05 critical value
t_crit = tinv(1 - alpha/2, n-2);

% calculate power - prob. of rejecting H0 = mass in H1 above t_crit
power = 100 * (1 - tcdf(t_crit - t1, n-2));

fprintf('With a sample size of %d, the power for detecting r >= 0.5 with a 2-tailed test is %.1f%%.\n', ...
    n, power);

% now find sample size needed for 99% power
n_99 = sampsizepwr('t', [0, se_r1 * sqrt(n)], r1, 0.99);

% check:
t2 = r1 / sqrt((1 - r1^2) / (n_99 - 2));
t_crit2 = tinv(1 - alpha/2, n_99-2);
power2 = 100 * (1 - tcdf(t_crit2 - t2, n_99-1));
assert(abs(power2 - 99) <= 1);

fprintf('A sample size of at least %d would provide a power of 99%%.\n', n_99);