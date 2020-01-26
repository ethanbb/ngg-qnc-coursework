%% Example 1: Quantal release
%% Exercise 1

n = 10; % 10 quanta available
p = 0.2; % release probability = 0.2
q = 1-p;

probK = @(k) nchoosek(n, k) * p^k * q^(n-k);

probs = zeros(11, 1);
for k = 0:10
    probs(k+1) = probK(k);
    fprintf('Pr(%d quanta): %.3f\n', k, probs(k+1));
end

figure;
bar(0:10, probs);
xlabel('k')
title('Probability of k quanta in 10 trials (p=0.2)');

%% Exercise 2

n = 14; % 14 available for release
k = 8;

ps = 0.1:0.1:1.0;

likelihood = binopdf(k, n, ps);
[prob, pk] = max(likelihood);

fprintf('The most likely p is %.1f with probability %.3f\n', ...
    ps(pk), prob);

%% Exercise 3

n = 14;
k1 = 8;
k2 = 5;

ps = 0.1:0.1:1.0;

likelihood1 = binopdf(k1, n, ps);
likelihood2 = binopdf(k2, n, ps);
likelihood = likelihood1 .* likelihood2;

loglike1 = log(likelihood1);
loglike2 = log(likelihood2);
loglike = loglike1 + loglike2;

fprintf('Likelihood of p=0.1: %.3f\n', likelihood(ps == 0.1));
fprintf('Log-likelihood of p=0.1: %.3f\n', loglike(ps == 0.1));

[~, kmax] = max(loglike);
fprintf('Most likely p at resolution of 0.1: %.1f\n', ps(kmax));

% Try resolution of 0.01
ps_hires = 0.01:0.01:1.0;
loglike1_hires = log(binopdf(k1, n, ps_hires));
loglike2_hires = log(binopdf(k2, n, ps_hires));
loglike_hires = loglike1_hires + loglike2_hires;

[~, kmax] = max(loglike_hires);
fprintf('Most likely p at resolution of 0.01: %.2f\n', ps_hires(kmax));

%% Exercise 4

% binofit expects scalar inputs and the documentation suggests using
% binofit(sum(k), n*length(k)) to fit to multiple samples, which
% makes sense since doing n independent trials m times is the same as
% just doing n*m independent trials.

n = 14;

ks = 0:14;
k_freqs = [0, 0, 3, 10, 19, 26, 16, 16, 5, 5, 0, 0, 0, 0, 0];

total_k = sum(ks .* k_freqs);
total_n = n * sum(k_freqs);
assert(total_n == n * 100, 'Oops there''s a mistake in k_freqs');

phat = binofit(total_k, total_n);
fprintf('Most likely p: %.2f\n', round(phat, 2));

%% Exercise 5

phat_new = binofit(7, 14); % (clearly this is just 0.5)
fprintf('New phat: %.2f\n', round(phat_new, 2));

p_if_null = binopdf(7, 14, 0.3);
fprintf('Prob. of result under null hypothesis (p) = %.4f\n', p_if_null);

% p = 0.0618 - not low enough to conclude there is an effect under any
% reasonable alpha. We would need to gather more data.

%% Bonus exercise

temps = [4.0, 3.5, 0.0, 2.0, 6.5, 3.0];
ntemps = length(temps);

ks = 0:4;
k_freqs = [
    615, 206, 33,  2,   0
    604, 339, 94,  11,  2
    332, 126, 21,  1,   0
    573, 443, 154, 28,  2
    172, 176, 89,  12,  1
    80,  224, 200, 32,  4
];

figure;

for kT = 1:ntemps
    subplot(3, 2, kT);
    hold on;
    
    fprintf('Temperature: %.1f C\n', temps(kT));
    quanta = repelem(ks, k_freqs(kT, :));
    v = var(quanta);
    m = mean(quanta);
    fprintf('Mean: %.2f quanta; variance: %.2f\n', m, v);
    p = 1 - v / m;
    fprintf('Estimated prob. of release: %.2f\n', p);
    n = m / p;
    fprintf('Estimated number of release sites: %.2f\n\n', n);

    plot(ks, k_freqs(kT, :) / sum(k_freqs(kT, :)), 'ro-');
    plot(ks, binopdf(ks, round(n), p), 'bo-');
    plot(ks, poisspdf(ks, m), 'go-');
    title(sprintf('Temp = %.1f, estimated n = %d', temps(kT), round(n)));
    legend('Empirical', 'Binomial fit', 'Poisson fit');
end

% Qualitatively, the empirical distributions do match the best-fitting
% binomial distributions very well. Interestingly, the best-fitting number
% of release sites depends on the temperature.