function compareModalities()
%COMPAREMODALITIES Quale sequenza RM usare: V, T2_FS o ADC (set dev, slice di area massima).
%   Per ogni paziente e sequenza: un seme al centro della maschera, ROI di 60 mm, normalizzazione e
%   filtro gaussiano, poi
%     CNR          contrasto tra tumore e anello di tessuto a 2-8 mm attorno (segno = polarita')
%     Dice Otsu    region growing dal seme con la soglia di Otsu
%     Dice oracolo miglior Dice su 46 soglie: limite superiore di una soglia unica
%   Uscite: results/compareModalities/per_patient.csv, boxplots.png
T = readCases('manifest/dev_cases.csv');
outDir = 'results/compareModalities'; [~, ~] = mkdir(outDir);
ROI_MM = 60; RING_MM = [2 8];
res = table();
for mo = unique(T.modality, 'stable')'
    for c = loadStartSlices(T, mo)'
        % ROI 60 x 60 mm centrata sul seme
        hx = round(ROI_MM / c.sp(1) / 2); hy = round(ROI_MM / c.sp(2) / 2);
        x0 = max(c.seed(1) - hx, 1); x1 = min(c.seed(1) + hx, size(c.img, 1));
        y0 = max(c.seed(2) - hy, 1); y1 = min(c.seed(2) + hy, size(c.img, 2));
        roi = c.img(x0:x1, y0:y1); ref = c.ref(x0:x1, y0:y1); s = c.seed - [x0 - 1, y0 - 1];
        % normalizzazione robusta (percentili 1-99) e filtro gaussiano (sigma 1 pixel)
        lo = prctile(roi(:), 1); hi = prctile(roi(:), 99);
        rn = imgaussfilt(min(max((roi - lo) / (hi - lo + eps), 0), 1), 1.0);
        % CNR tumore contro anello di 2-8 mm attorno alla maschera
        px = @(mm) max(round(mm / c.sp(1)), 1);
        ring = imdilate(ref, strel('disk', px(RING_MM(2)))) & ~imdilate(ref, strel('disk', px(RING_MM(1))));
        cnr = (mean(rn(ref)) - mean(rn(ring))) / (std(rn(ring)) + eps);
        % Dice con Otsu e con la soglia oracolo, per tumore chiaro e scuro
        t = graythresh(rn); dOtsu = [0 0]; dOracle = [0 0]; pols = ["bright", "dark"];
        for p = 1:2
            dOtsu(p) = diceCoef(regionGrow(rn, s, t, pols(p), px(3)), ref);
            for q = linspace(0.05, 0.95, 46)
                dOracle(p) = max(dOracle(p), diceCoef(regionGrow(rn, s, q, pols(p), px(3)), ref));
            end
        end
        res = [res; table(c.id, mo, c.k, cnr, t, dOtsu(1), dOtsu(2), dOracle(1), dOracle(2), 'VariableNames', ...
            {'patient', 'modality', 'slice', 'cnr', 'otsu_t', 'dice_otsu_bright', 'dice_otsu_dark', ...
             'dice_oracle_bright', 'dice_oracle_dark'})]; %#ok<AGROW>
    end
end
writetable(res, fullfile(outDir, 'per_patient.csv'));

%% Sintesi per sequenza: polarita' = segno del CNR mediano
for mo = unique(res.modality, 'stable')'
    x = res(res.modality == mo, :); pol = "bright"; if median(x.cnr) < 0, pol = "dark"; end
    fprintf('%-6s polarita'' %-6s CNR med %.2f | Dice Otsu med %.3f | Dice oracolo med %.3f\n', mo, pol, ...
        median(x.cnr), median(x.("dice_otsu_" + pol)), median(x.("dice_oracle_" + pol)));
end
f = figure('Visible', 'off'); theme(f, 'light');
boxchart(categorical(res.modality), res.cnr); yline(0, ':'); ylabel('CNR (tumour vs 2-8 mm ring)');
exportgraphics(f, fullfile(outDir, 'boxplots.png')); close(f);
end

function bw = regionGrow(rn, s, t, pol, r)
% Region growing = ricostruzione morfologica del seme dentro la maschera sogliata
% (equivale alla componente connessa che contiene il seme). Se il pixel del seme e' fuori soglia
% si parte dai pixel sopra soglia entro r pixel, come seedFallback nel metodo.
    if pol == "bright", M = rn > t; else, M = rn < t; end
    seedIm = false(size(M)); seedIm(s(1), s(2)) = true;
    if ~M(s(1), s(2)), seedIm = imdilate(seedIm, strel('disk', r, 0)) & M; end
    bw = imreconstruct(seedIm, M);
end
