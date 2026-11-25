function [p, z, n] = wilcoxonSignedRank(d)
%WILCOXONSIGNEDRANK Test di Wilcoxon per dati appaiati (differenze d), approssimazione normale con
%   correzione per i pareggi. Scritto a mano perche' signrank richiede lo Statistics Toolbox.
    d = d(:); d = d(~isnan(d) & d ~= 0); n = numel(d);
    if n < 5, p = NaN; z = NaN; return; end
    a = abs(d); [s, idx] = sort(a); r = zeros(n, 1); t = []; i = 1;
    while i <= n
        j = i; while j < n && s(j + 1) == s(i), j = j + 1; end
        r(idx(i:j)) = (i + j) / 2; if j > i, t(end+1) = j - i + 1; end %#ok<AGROW>
        i = j + 1;
    end
    W = sum(r(d > 0)); mu = n * (n + 1) / 4;
    sigma = sqrt(n * (n + 1) * (2 * n + 1) / 24 - sum(t.^3 - t) / 48);
    z = (W - mu) / sigma; p = erfc(abs(z) / sqrt(2));
end
