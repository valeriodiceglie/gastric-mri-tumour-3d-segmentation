function seed = centroidSeed(B)
%CENTROIDSEED Pixel della maschera B piu' vicino al suo baricentro, come [riga colonna].
%   Simula il clic dell'operatore al centro di una regione. Si prende un pixel della maschera,
%   e non il baricentro stesso, perche' in una regione a forma di C il baricentro cade fuori.
    [xs, ys] = find(B);
    [~, j] = min((xs - mean(xs)).^2 + (ys - mean(ys)).^2);
    seed = [xs(j) ys(j)];
end
