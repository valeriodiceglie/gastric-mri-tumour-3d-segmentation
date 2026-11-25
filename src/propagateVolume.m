function M3 = propagateVolume(I, sp, k0, seeds0, P, Q)
%PROPAGATEVOLUME Segmenta la slice di partenza k0 e propaga la maschera alle slice vicine.
%   Verso l'alto e poi verso il basso: la maschera della slice precedente fa da seme sulla
%   successiva (segmentSlice), finche' un criterio di arresto non scatta.
%
%   I      : volume [X Y Z];  sp : spaziatura [sx sy sz] (mm)
%   k0     : slice di partenza;  seeds0 : semi su k0 (N x 2 [riga colonna])
%   P      : parametri di segmentSlice
%   Q      : parametri della propagazione (valori finali in methodConfig)
%     marker          'points' = un seme per pezzo della maschera precedente (pickSeeds)
%                     'mask'   = la maschera precedente (erosa) e' il seme e il campione di intensita'
%     erodeMM         erosione della maschera precedente prima di usarla come seme (mm)
%     clipMM          la nuova maschera resta entro la precedente dilatata di clipMM (0 = nessun vincolo)
%     minSeedAreaMM2, maxSeeds   regola dei semi in modalita' 'points' (pickSeeds)
%     arresto quando la nuova maschera:
%     minAreaMM2      ha area sotto questa soglia (mm^2)
%     minOverlap      ha Dice con la precedente sotto questa soglia (0 = criterio spento)
%   Campi opzionali (se mancano, il criterio e' spento):
%     maxGrowth       area nuova <= maxGrowth x area precedente; se la supera, il vincolo attorno
%                     alla maschera precedente si stringe finche' rientra (Inf = spento)
%     roiFollow       true = la ROI copre la maschera precedente + roiMarginMM per lato
%     minRelArea      arresto: area sotto questa frazione dell'area massima raggiunta (0 = spento)
%     minThickMM      arresto: spessore sotto questa soglia (wallThickness; 0 = spento)
%     thickStat       statistica dello spessore: 'max' | 'skel'
%   M3     : maschera 3D logica
    Q = withDefaults(Q);
    pxmm = mean(sp(1:2)); nz = size(I, 3);
    M3 = false(size(I));
    M0 = segmentSlice(double(I(:, :, k0)), seeds0, sp(1:2), P);
    M3(:, :, k0) = M0;
    if ~any(M0(:)), return; end
    Pk = P; Pk.perSeed = false; Pk.roiRef = [];
    if Q.roiFollow, Pk.roiMM = 2 * Q.roiMarginMM; end
    aMax = nnz(M0);
    for dir = [1 -1]
        prev = M0; k = k0 + dir;
        while k >= 1 && k <= nz
            % ---- seme dalla maschera precedente ----
            mk = prev; re = round(Q.erodeMM / pxmm);
            if re >= 1
                e = imerode(prev, strel('disk', re, 0));
                if any(e(:)), mk = e; end              % se l'erosione cancella tutto, si usa la maschera intera
            end
            if strcmp(Q.marker, 'points')
                sd = pickSeeds(mk, sp(1:2), Q.minSeedAreaMM2);
                sd = sd(1:min(end, Q.maxSeeds), :);
            else
                sd = mk;
            end
            % ---- segmentazione vincolata attorno alla maschera precedente ----
            if Q.roiFollow, Pk.roiRef = prev; end
            m = segmentSlice(double(I(:, :, k)), sd, sp(1:2), Pk);
            rc = round(Q.clipMM / pxmm);
            if rc >= 1, m = m & imdilate(prev, strel('disk', rc, 0)); end
            aPrev = nnz(prev);
            if nnz(m) > Q.maxGrowth * aPrev            % limite di crescita: vincolo via via piu' stretto
                ok = false;
                for r = max(rc, round(4 / pxmm)) - 1 : -1 : 1
                    m2 = m & imdilate(prev, strel('disk', r, 0));
                    if nnz(m2) <= Q.maxGrowth * aPrev, m = m2; ok = true; break; end
                end
                if ~ok, m = m & prev; end
            end
            % ---- criteri di arresto ----
            if nnz(m) * sp(1) * sp(2) < Q.minAreaMM2, break; end
            if nnz(m) < Q.minRelArea * aMax, break; end
            if Q.minThickMM > 0 && wallThickness(m, pxmm, Q.thickStat) < Q.minThickMM, break; end
            if Q.minOverlap > 0 && diceCoef(m, prev) < Q.minOverlap, break; end
            M3(:, :, k) = m; prev = m; aMax = max(aMax, nnz(m));
            k = k + dir;
        end
    end
end

function Q = withDefaults(Q)
% Campi opzionali mancanti: criterio spento.
    d = struct('maxGrowth', Inf, 'roiFollow', false, 'roiMarginMM', 20, 'minRelArea', 0, ...
               'minThickMM', 0, 'thickStat', 'max');
    for f = fieldnames(d)'
        if ~isfield(Q, f{1}), Q.(f{1}) = d.(f{1}); end
    end
end
