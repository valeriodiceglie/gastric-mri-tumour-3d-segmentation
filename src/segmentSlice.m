function [mask, box, steps] = segmentSlice(I2, seeds, sp, P)
%SEGMENTSLICE Segmenta il tumore su una slice a partire dai semi (region growing).
%   Passi: ROI attorno ai semi -> normalizzazione (percentili 1-99) -> filtro gaussiano ->
%   soglia (Otsu, o adattiva attorno ai semi) -> region growing dai semi -> apertura ->
%   chiusura e riempimento dei buchi.
%
%   I2    : slice 2D (double)
%   seeds : N x 2 [riga colonna] dei semi, dentro il tumore; oppure una maschera logica grande
%           come I2 (propagazione: la regione fa da seme e da campione di intensita')
%   sp    : spaziatura [sx sy] (mm)
%   P     : parametri (valori finali in methodConfig)
%     roiMM          lato della ROI quadrata (mm); con piu' semi copre tutti i semi + roiMM/2 per lato
%     sigmaMM        sigma del filtro gaussiano (mm; 0 = nessun filtro)
%     method         'otsu' (soglia globale nella ROI) | 'adaptive' (media e dev. std. vicino ai semi)
%     sided, k       solo 'adaptive': 'one' = pixel non piu' scuri di media - k dev. std. (tumore
%                    chiaro; il contrario se scuro), 'two' = pixel entro media +- k dev. std.
%     seedRadiusMM   raggio dell'intorno dei semi da cui si stimano media e dev. std. (mm)
%     polarity       +1 tumore chiaro (V), -1 tumore scuro (ADC)
%     openMM         raggio dell'apertura (mm; 0 = non applicata)
%     closeMM        raggio della chiusura (mm; 0 = non applicata)
%     fill           riempimento dei buchi (true/false)
%     conn           connettivita' del region growing (4 o 8)
%     perSeed        (opzionale) true = soglia e crescita separate per ogni seme, poi unione
%     seedFallback   (opzionale) true = un seme che cade su un pixel sotto soglia e' sostituito dai
%                    pixel sopra soglia entro seedRadiusMM (altrimenti la crescita non parte)
%     roiRef         (opzionale) maschera logica: la ROI copre anche il suo riquadro
%   mask  : maschera logica grande come I2
%   box   : riquadro della ROI [x0 x1 y0 y1] (per ritagliare le figure)
%   steps : immagini intermedie nella ROI, per le figure: roi (normalizzata e filtrata),
%           threshold (maschera sogliata), grown (region growing prima dell'apertura)
%           (vuote con perSeed)
    [nx, ny] = size(I2);
    pxmm = mean(sp(1:2));
    maskMode = islogical(seeds) && isequal(size(seeds), size(I2));
    if maskMode, [px, py] = find(seeds); pts = [px py]; else, pts = seeds; end

    % ---- 1. ROI quadrata attorno ai semi (e alla maschera P.roiRef, se data) ----
    boxPts = pts;
    if isfield(P, 'roiRef') && ~isempty(P.roiRef) && any(P.roiRef(:))
        [rx, ry] = find(P.roiRef); boxPts = [rx ry; pts];
    end
    hx = round(P.roiMM / sp(1) / 2);  hy = round(P.roiMM / sp(2) / 2);
    x0 = max(min(boxPts(:, 1)) - hx, 1); x1 = min(max(boxPts(:, 1)) + hx, nx);
    y0 = max(min(boxPts(:, 2)) - hy, 1); y1 = min(max(boxPts(:, 2)) + hy, ny);
    roi = I2(x0:x1, y0:y1);
    s = pts - [x0 - 1, y0 - 1];                              % semi nelle coordinate della ROI

    % ---- 2. normalizzazione robusta e filtro gaussiano ----
    lo = prctile(roi(:), 1); hi = prctile(roi(:), 99);
    rn = min(max((roi - lo) / (hi - lo + eps), 0), 1);
    if P.sigmaMM > 0, rn = imgaussfilt(rn, P.sigmaMM / pxmm); end

    % ---- 3. soglia, region growing e apertura: tutti i semi insieme, oppure seme per seme ----
    perSeed = ~maskMode && isfield(P, 'perSeed') && P.perSeed && size(s, 1) > 1;
    if maskMode
        seedIm = seeds(x0:x1, y0:y1);
        M = thresholdMask(rn, seedIm, P, pxmm, true);
        [G, grown] = growAndOpen(seedIm, M, P, pxmm);
    elseif ~perSeed
        seedIm = false(size(rn)); seedIm(sub2ind(size(rn), s(:, 1), s(:, 2))) = true;
        M = thresholdMask(rn, seedIm, P, pxmm, false);
        [G, grown] = growAndOpen(fallbackSeeds(seedIm, M, P, pxmm), M, P, pxmm);
    else
        G = false(size(rn)); M = []; grown = [];
        for q = 1:size(s, 1)
            sq = false(size(rn)); sq(s(q, 1), s(q, 2)) = true;
            Mq = thresholdMask(rn, sq, P, pxmm, false);
            G = G | growAndOpen(fallbackSeeds(sq, Mq, P, pxmm), Mq, P, pxmm);
        end
    end

    % ---- 4. chiusura e riempimento dei buchi sul risultato complessivo ----
    rc = round(P.closeMM / pxmm);
    if rc >= 1, G = imclose(G, strel('disk', rc, 0)); end
    if P.fill, G = imfill(G, 'holes'); end
    mask = false(nx, ny); mask(x0:x1, y0:y1) = G;
    box = [x0 x1 y0 y1];
    steps = struct('roi', rn, 'threshold', M, 'grown', grown);
end

function M = thresholdMask(rn, seedIm, P, pxmm, useRegion)
% Soglia: Otsu sulla ROI, oppure intervallo attorno alla media del tessuto vicino ai semi.
% useRegion = true: il campione e' la regione-seme stessa (propagazione), non un disco attorno ai semi.
    switch P.method
        case 'otsu'
            t = graythresh(rn);
            if P.polarity > 0, M = rn > t; else, M = rn < t; end
        case 'adaptive'
            if useRegion
                nb = seedIm;
            else
                r  = max(round(P.seedRadiusMM / pxmm), 1);
                nb = imdilate(seedIm, strel('disk', r, 0));
            end
            mu = mean(rn(nb)); sd = max(std(rn(nb)), 0.01);
            if strcmp(P.sided, 'two')
                M = abs(rn - mu) <= P.k * sd;
            elseif P.polarity > 0
                M = rn >= mu - P.k * sd;          % tumore chiaro: tutto cio' che non e' piu' scuro
            else
                M = rn <= mu + P.k * sd;          % tumore scuro: tutto cio' che non e' piu' chiaro
            end
    end
end

function s = fallbackSeeds(s, M, P, pxmm)
% Un seme che cade su un pixel piu' scuro della soglia (rumore, zona necrotica) non fa crescere
% nulla e il paziente resta vuoto. Con P.seedFallback quel seme e' sostituito dai pixel sopra soglia
% entro seedRadiusMM; i semi sopra soglia restano come sono.
    if ~(isfield(P, 'seedFallback') && P.seedFallback), return; end
    dark = s & ~M;
    if any(dark(:))
        r = max(round(P.seedRadiusMM / pxmm), 1);
        s = (s & M) | (imdilate(dark, strel('disk', r, 0)) & M);
    end
end

function [G, grown] = growAndOpen(seedIm, M, P, pxmm)
% Region growing = ricostruzione morfologica dei semi dentro la maschera sogliata; poi apertura.
% grown = regione prima dell'apertura.
    G = imreconstruct(seedIm & M, M, P.conn); grown = G;
    ro = round(P.openMM / pxmm);
    if ro >= 1 && any(G(:))
        se = strel('disk', ro, 0);
        G2 = imopen(G, se);
        mk = imdilate(seedIm, se) & G2;          % riparte dai semi: taglia i "ponti" verso altri tessuti
        if any(mk(:)), G = imreconstruct(mk, G2, P.conn); end   % se l'apertura cancella tutto, si tiene G
    end
end
