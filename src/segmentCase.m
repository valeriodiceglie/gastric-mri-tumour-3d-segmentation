function [M3, k0, seeds] = segmentCase(I, S, sp, P, Q)
%SEGMENTCASE Metodo completo su un paziente (modalita' V, parametri di methodConfig).
%   1) slice di partenza = slice di area massima della maschera esperta S;
%   2) un seme per ogni pezzo del tumore su quella slice (simula i clic dell'operatore, pickSeeds);
%   3) segmentazione della slice e propagazione alle vicine (propagateVolume).
%
%   I  : volume [X Y Z];  S : maschera esperta logica;  sp : spaziatura [sx sy sz] (mm)
%   P, Q (opzionali): parametri al posto di quelli di methodConfig (varianti, es. senza seedFallback)
%   M3 : maschera 3D;  k0 : slice di partenza;  seeds : semi su k0 (N x 2)
    if nargin < 4, [P, Q] = methodConfig(); end
    [~, k0] = max(squeeze(sum(S, [1 2])));
    seeds = pickSeeds(S(:, :, k0), sp(1:2), Q.minSeedAreaMM2);
    M3 = propagateVolume(I, sp, k0, seeds, P, Q);
end
