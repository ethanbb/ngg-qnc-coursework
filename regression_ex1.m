%% Exercise 1: regression on age/wing length

age = [3; 4; 5; 6; 7; 8; 9; 11; 12; 14; 15; 16; 17];
wlen = [1.4; 1.5; 2.2; 2.4; 3.1; 3.2; 3.2; 3.9; 4.1; 4.7; 4.5; 5.2; 5];

n = length(age);
assert(n == length(wlen));

%% 1. Plot

figure;

plot(age, wlen, '*');
xlabel('Age');
ylabel('Wing length');

%% 2. Find regression line

X = [ones(n, 1), age];
% equation: wlen = X * [a; b]
ab = X \ wlen;
a = ab(1); b = ab(2);

hold on;
h1 = plot(age, a + b * age, 'k--');
legend(h1, sprintf('Linear fit: y = %f + %f * age', a, b), 'Location', 'southeast');

%% 3. Can we reject H0: b = 0?

% It sure looks like we can, but find the SE and do a t-test

residuals = wlen - (a + b * age);
seb = norm(residuals) / norm(age - mean(age)) / sqrt(n-2);
tstat = b / seb;
p = 2 * (1 - tcdf(abs(tstat), n - 2));

fprintf('p-value for b != 0 is %g.\n', p);
alpha = 0.05;
if p < alpha
    disp('This is significant; we reject H0.');
else
    disp('We cannot reject H0.');
end

%% 4. Confidence intervals

t_alpha = tinv([alpha/2, 1-alpha/2], n-2);
ci = b + t_alpha * seb;
fprintf('Confidence interval for b: [%f, %f]\n', ci(1), ci(2));

h2 = plot(age, a + ci(1) * age, 'k:');
plot(age, a + ci(2) * age, 'k:');

legend([h1, h2], get(gca, 'Legend').String{1}, '95% CI');

%% 5. To get r, multiply by std dev of X and divide by std dev of Y:

r = b * std(age) / std(wlen);
fprintf('Age and wing length have an r of %.2f.\n', r);

% in order for the t test to work the same way, the std error
% of r must also be scaled in this way:
ser = norm(residuals) / norm(wlen - mean(wlen)) / sqrt(n-2);
assert(r / ser == b / seb);

%% 6. Try adding some noise

wlen2 = wlen + 3*randn(size(wlen));
ab2 = X \ wlen2;
a2 = ab2(1); b2 = ab2(2);

plot(age, wlen2, 'r*');
plot(age, a2 + b2 * age, 'r--');

% b2 is still close to b, but offset.
resid2 = wlen2 - (a2 + b2 * age);
seb2 = norm(resid2) / norm(age - mean(age)) / sqrt(n-2);
tstat2 = b2 / seb2;

p2 = 2 * (1 - tcdf(abs(tstat2), n - 2));

% p is now 0.38 - definitely not significant.