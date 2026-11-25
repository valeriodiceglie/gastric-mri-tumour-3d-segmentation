function C = loadVolumes(T, modality, zMargin)
%LOADVOLUMES Volumi di una modalita' in memoria, con slice di partenza e semi.
%   T        : tabella dei casi (readCases)
%   modality : "V" o "ADC"
%   zMargin  : (opzionale) tiene solo le slice del tumore +- zMargin (default Inf = volume intero)
%   C(i)     : id, vol (int16, per contenere la memoria), ref (maschera esperta), sp [sx sy sz] (mm),
%              k0 (slice di area massima), seeds (un seme per pezzo su k0, pickSeeds),
%              nTumorSlices (slice con tumore)
%   Pazienti in ordine di ID.
    if nargin < 3, zMargin = Inf; end
    [~, Q] = methodConfig();
    T = sortrows(T(T.modality == modality, :), 'patient_id');
    for i = height(T):-1:1
        [I, S, sp] = readCase(T.img(i), T.seg(i));
        area = squeeze(sum(S, [1 2])); ks = find(area);
        [~, k0] = max(area);
        z = max(ks(1) - zMargin, 1):min(ks(end) + zMargin, size(I, 3));
        C(i) = struct('id', T.patient_id(i), 'vol', int16(I(:, :, z)), 'ref', S(:, :, z), 'sp', sp, ...
                      'k0', k0 - z(1) + 1, 'seeds', pickSeeds(S(:, :, k0), sp(1:2), Q.minSeedAreaMM2), ...
                      'nTumorSlices', numel(ks));
    end
    C = C(:);
end
