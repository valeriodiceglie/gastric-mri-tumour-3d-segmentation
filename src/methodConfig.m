function [P, Q, Qv1] = methodConfig(modality)
%METHODCONFIG Parametri del metodo, scelti sul set dev e poi tenuti fissi sul test.
%   [P, Q] = methodConfig()  o  methodConfig("V") : metodo finale (modalita' V)
%     P   : segmentazione di una slice (segmentSlice)
%     Q   : propagazione tra slice (propagateVolume)
%     Qv1 : prima versione della propagazione (senza limite di crescita ne' arresto relativo),
%           tenuta solo come confronto nelle ablazioni e come base della sua griglia
%   P = methodConfig("ADC") : segmentazione della slice di partenza scelta per ADC (solo confronto)
%   seedFallback: aggiunto dopo un primo passaggio sul test (16 risultati vuoti; lo stesso difetto in
%   1 paziente dev). Non ha parametri; il primo passaggio e' confrontato in analyzeSeedFallback.
%   Non modificare guardando il test: ogni modifica va valutata sul dev e dichiarata come post hoc.
    if nargin < 1, modality = "V"; end
    switch modality
        case "V"
            P = struct('method', 'adaptive', 'sided', 'one', 'k', 1, 'roiMM', 60, 'sigmaMM', 0.5, 'openMM', 3, ...
                       'closeMM', 0, 'fill', true, 'conn', 4, 'polarity', 1, 'seedRadiusMM', 3, 'perSeed', false, ...
                       'seedFallback', true);
        case "ADC"
            P = struct('method', 'adaptive', 'sided', 'two', 'k', 2, 'roiMM', 60, 'sigmaMM', 1.5, 'openMM', 0, ...
                       'closeMM', 2, 'fill', true, 'conn', 4, 'polarity', -1, 'seedRadiusMM', 3, 'perSeed', false, ...
                       'seedFallback', true);
    end
    Q = struct('marker', 'mask', 'erodeMM', 1, 'clipMM', 1, 'minAreaMM2', 20, 'minOverlap', 0.2, 'minSeedAreaMM2', 10, ...
               'maxSeeds', Inf, 'maxGrowth', 1.2, 'roiFollow', false, 'roiMarginMM', 20, 'minRelArea', 0.15);
    Qv1 = struct('marker', 'points', 'erodeMM', 1, 'clipMM', 2, 'minAreaMM2', 20, 'minOverlap', 0.2, ...
                 'minSeedAreaMM2', 10, 'maxSeeds', Inf);
end
