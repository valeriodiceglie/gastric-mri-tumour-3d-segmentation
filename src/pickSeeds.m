function seeds = pickSeeds(B, sp, minAreaMM2)
%PICKSEEDS Un seme per ogni pezzo (componente connessa) della maschera B.
%   Simula un operatore che clicca al centro di ogni pezzo visibile del tumore: la maschera
%   esperta e' spesso una parete gastrica divisa in piu' pezzi, e un solo seme ne raggiunge uno.
%
%   B          : maschera 2D logica
%   sp         : spaziatura [sx sy] (mm)
%   minAreaMM2 : i pezzi piu' piccoli sono scartati (il piu' grande e' sempre tenuto)
%   seeds      : N x 2 [riga colonna], dal pezzo piu' grande al piu' piccolo (centroidSeed)
    L = bwlabel(B, 8);
    areas = zeros(max(L(:)), 1);
    for c = 1:numel(areas), areas(c) = nnz(L == c) * sp(1) * sp(2); end
    [~, order] = sort(areas, 'descend');
    seeds = zeros(0, 2);
    for c = order'
        if areas(c) < minAreaMM2 && c ~= order(1), continue; end
        seeds(end+1, :) = centroidSeed(L == c); %#ok<AGROW>
    end
end
