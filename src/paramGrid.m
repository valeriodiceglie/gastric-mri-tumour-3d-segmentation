function G = paramGrid(base, varargin)
%PARAMGRID Griglia di parametri: tutte le combinazioni dei valori dati, come array di struct.
%   G = paramGrid(base, 'k', [1 2], 'method', {'otsu', 'adaptive'})
%   parte da base e assegna ai campi elencati ogni combinazione di valori (testo in una cella).
%   L'ordine e' quello di cicli for annidati: il primo campo e' il ciclo esterno, l'ultimo varia
%   piu' in fretta. A parita' di Dice vince la prima configurazione, quindi l'ordine conta.
    names = varargin(1:2:end); vals = varargin(2:2:end);
    n = cellfun(@numel, vals);
    G = repmat(base, 1, prod(n));
    for g = 1:prod(n)
        sub = cell(1, numel(n));
        [sub{end:-1:1}] = ind2sub(fliplr(n), g);   % l'ultimo campo e' l'indice che varia piu' in fretta
        for f = 1:numel(names)
            v = vals{f};
            if iscell(v), G(g).(names{f}) = v{sub{f}}; else, G(g).(names{f}) = v(sub{f}); end
        end
    end
end
