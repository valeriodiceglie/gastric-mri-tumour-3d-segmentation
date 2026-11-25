function C = loadStartSlices(T, modality)
%LOADSTARTSLICES Slice di partenza (area massima della maschera esperta) di ogni paziente.
%   T        : tabella dei casi (readCases)
%   modality : "V", "T2_FS" o "ADC"
%   C(i)     : id, img (slice 2D), ref (maschera esperta sulla slice), sp [sx sy] (mm),
%              k (indice della slice), seed (un seme al centro della maschera, centroidSeed)
%   Pazienti in ordine di ID.
    T = sortrows(T(T.modality == modality, :), 'patient_id');
    for i = height(T):-1:1
        [I, S, sp] = readCase(T.img(i), T.seg(i));
        [~, k] = max(squeeze(sum(S, [1 2])));
        C(i) = struct('id', T.patient_id(i), 'img', I(:, :, k), 'ref', S(:, :, k), 'sp', sp(1:2), ...
                      'k', k, 'seed', centroidSeed(S(:, :, k)));
    end
    C = C(:);
end
