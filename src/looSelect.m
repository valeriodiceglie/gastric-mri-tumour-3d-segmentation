function [loo, pick] = looSelect(D)
%LOOSELECT Leave-one-out della taratura. D : configurazioni x pazienti (Dice).
%   Per ogni paziente i sceglie la configurazione con la mediana migliore sugli altri pazienti
%   (pick(i)) e restituisce il suo Dice su i (loo(i)): stima onesta del risultato della taratura.
    n = size(D, 2); loo = zeros(n, 1); pick = zeros(n, 1);
    for i = 1:n
        [~, pick(i)] = max(median(D(:, setdiff(1:n, i)), 2));
        loo(i) = D(pick(i), i);
    end
end
