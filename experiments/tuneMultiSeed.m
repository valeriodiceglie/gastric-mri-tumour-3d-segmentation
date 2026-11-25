function tuneMultiSeed()
%TUNEMULTISEED Piu' semi sulla slice di partenza: uno per ogni pezzo del tumore (set dev, V e ADC).
%   La maschera esperta e' spesso una parete gastrica in piu' pezzi: un seme solo ne raggiunge uno.
%   Confronta, con Dice sull'intera slice:
%     1) un seme (configurazione di tuneStartSlice)
%     2) un seme per pezzo (pickSeeds), stessa configurazione, crescita unica da tutti i semi
%     3) un seme per pezzo, stessa configurazione, soglia e crescita separate per ogni seme
%     4) un seme per pezzo, parametri ritarati (griglia di 1152 configurazioni attorno alla 1)
%     5) leave-one-out della ritaratura (stima onesta)
%   Uscite in results/tuneMultiSeed/: grid_, summary_, per_patient_<mod>.csv, contact_best_<mod>.png
T = readCases('manifest/dev_cases.csv');
outDir = 'results/tuneMultiSeed'; [~, ~] = mkdir(outDir); tic;
MIN_AREA = [10 30];                       % area minima (mm^2) di un pezzo per ricevere un seme
names = {'k', 'roiMM', 'sigmaMM', 'openMM', 'closeMM', 'fill', 'perSeed', 'minAreaMM2'};

for mo = ["V", "ADC"]
    B = methodConfig(mo); B.minAreaMM2 = MIN_AREA(1);     % configurazione scelta con un seme
    C = loadStartSlices(T, mo); n = numel(C);
    seeds = cell(n, numel(MIN_AREA));
    for i = 1:n
        for a = 1:numel(MIN_AREA), seeds{i, a} = pickSeeds(C(i).ref, C(i).sp, MIN_AREA(a)); end
    end
    nSeeds = cellfun(@(s) size(s, 1), seeds);
    diceOf = @(i, sd, P) diceCoef(segmentSlice(C(i).img, sd, C(i).sp, P), C(i).ref);

    %% Confronti a configurazione fissa
    Bper = B; Bper.perSeed = true;
    d1 = zeros(n, 1); d2 = d1; d3 = d1;
    for i = 1:n
        d1(i) = diceOf(i, C(i).seed, B);
        d2(i) = diceOf(i, seeds{i, 1}, B);
        d3(i) = diceOf(i, seeds{i, 1}, Bper);
    end

    %% Griglia attorno alla configurazione scelta con un seme
    cfg = paramGrid(B, 'k', B.k * [0.75 1 1.25], 'sigmaMM', [0.5 1 1.5], 'openMM', [0 1 2 3], ...
        'closeMM', [0 2], 'fill', [false true], 'roiMM', [60 80], 'perSeed', [false true], 'minAreaMM2', MIN_AREA);
    D = zeros(numel(cfg), n);
    for g = 1:numel(cfg)
        a = find(MIN_AREA == cfg(g).minAreaMM2);
        for i = 1:n, D(g, i) = diceOf(i, seeds{i, a}, cfg(g)); end
        if mod(g, 300) == 0, fprintf('%s: %d/%d configurazioni (%.0f s)\n', mo, g, numel(cfg), toc); end
    end
    med = median(D, 2); [~, best] = max(med);
    [loo, looCfg] = looSelect(D);
    R = struct2table(cfg); R.dice_median = med; R.dice_mean = mean(D, 2);
    R.dice_q25 = prctile(D, 25, 2); R.dice_q75 = prctile(D, 75, 2);
    writetable(sortrows(R, 'dice_median', 'descend'), fullfile(outDir, "grid_" + mo + ".csv"));

    %% Riepilogo (dice_median_multi = solo i pazienti con piu' pezzi)
    variants = ["Un seme"; "Piu' semi, crescita unica"; "Piu' semi, crescita per seme"; ...
                "Piu' semi, config ritarata (migliore)"; "Leave-one-out della ritaratura"];
    X = [d1'; d2'; d3'; D(best, :); loo'];
    cfgTxt = [describeParams(B, names); describeParams(B, names); describeParams(Bper, names); ...
              describeParams(cfg(best), names); sprintf('%d config diverse nei %d fold', numel(unique(looCfg)), n)];
    multi = nSeeds(:, 1) > 1;
    Sm = table(variants, cfgTxt, median(X, 2), mean(X, 2), prctile(X, 25, 2), prctile(X, 75, 2), median(X(:, multi), 2), ...
        'VariableNames', {'variante', 'configurazione', 'dice_median', 'dice_mean', 'dice_q25', 'dice_q75', 'dice_median_multi'});
    writetable(Sm, fullfile(outDir, "summary_" + mo + ".csv"));
    writetable(table([C.id]', nSeeds(:, 1), nSeeds(:, 2), d1, d2, d3, D(best, :)', loo, 'VariableNames', ...
        {'patient_id', 'n_semi_min10', 'n_semi_min30', 'dice_un_seme', 'dice_multi_unica', 'dice_multi_perseme', ...
         'dice_best', 'dice_loo'}), fullfile(outDir, "per_patient_" + mo + ".csv"));
    fprintf('\n== %s ==  pazienti con piu'' pezzi (>= 10 mm^2): %d su %d; semi per paziente: mediana %g, max %d\n', ...
        mo, nnz(multi), n, median(nSeeds(:, 1)), max(nSeeds(:, 1)));
    disp(Sm);

    %% Figura: configurazione migliore su ogni paziente
    P = cfg(best); a = find(MIN_AREA == P.minAreaMM2);
    for i = n:-1:1
        sd = seeds{i, a}; [m, box] = segmentSlice(C(i).img, sd, C(i).sp, P); id = char(C(i).id);
        tiles(i) = struct('img', C(i).img, 'ref', C(i).ref, 'mask', m, 'box', box, 'seeds', sd, ...
                          'label', sprintf('%s  %d seeds  Dice %.2f', id(end-5:end), size(sd, 1), D(best, i)));
    end
    saveContactSheet(fullfile(outDir, "contact_best_" + mo + ".png"), tiles, mo + " - " + describeParams(P, names));
    clear tiles
end
fprintf('\nFatto in %.0f s. Risultati in %s\n', toc, outDir);
end
