function ci = bootMedianCI(x, B, seed)
%BOOTMEDIANCI Intervallo di confidenza al 95% della mediana con bootstrap percentile sui pazienti.
    if nargin < 2, B = 2000; end; if nargin < 3, seed = 42; end
    x = x(~isnan(x)); n = numel(x); rs = RandStream('mt19937ar', 'Seed', seed);
    m = zeros(B, 1); for b = 1:B, m(b) = median(x(randi(rs, n, n, 1))); end
    ci = prctile(m, [2.5 97.5]);
end
