function tuneStartSlice()
%TUNESTARTSLICE Taratura della segmentazione della slice di partenza, con un seme (set dev, V e ADC).
%   Griglia di 2112 configurazioni di segmentSlice: soglia di Otsu o adattiva attorno al seme,
%   ROI, filtro gaussiano e morfologia. Dice sull'intera slice (non solo nella ROI).
%   Stima onesta della taratura con leave-one-out: per ogni paziente la configurazione e' scelta
%   sugli altri 19 e poi applicata a lui (looSelect).
%   Uscite in results/tuneStartSlice/, per modalita':
%     grid_<mod>.csv         tutte le configurazioni, ordinate per Dice mediano
%     summary_<mod>.csv      Otsu di partenza, Otsu tarato, miglior configurazione, leave-one-out
%     per_patient_<mod>.csv  Dice per paziente delle stesse quattro varianti
%     contact_best_<mod>.png risultato della miglior configurazione su ogni paziente
T = readCases('manifest/dev_cases.csv');
outDir = 'results/tuneStartSlice'; [~, ~] = mkdir(outDir); tic;
polarity = struct('V', 1, 'ADC', -1);

%% Griglia: Otsu (senza k) e soglia adattiva (da un lato o dai due lati, k deviazioni standard)
base = struct('method', '', 'sided', '', 'k', NaN, 'roiMM', 0, 'sigmaMM', 0, 'openMM', 0, ...
              'closeMM', 0, 'fill', false, 'conn', 0, 'seedFallback', true);
common = {'roiMM', [60 80], 'sigmaMM', [0.5 1 1.5], 'openMM', [0 1 2 3], 'closeMM', [0 2], ...
          'fill', [false true], 'conn', [4 8]};
cfg = [paramGrid(base, 'method', {'otsu'}, 'sided', {'-'}, 'k', NaN, common{:}), ...
       paramGrid(base, 'method', {'adaptive'}, 'sided', {'one', 'two'}, 'k', [1 1.5 2 2.5 3], common{:})];
names = ['method', 'sided', 'k', common(1:2:end)];
isOtsu = strcmp({cfg.method}, 'otsu');
otsuStart = find(isOtsu & [cfg.roiMM] == 60 & [cfg.sigmaMM] == 1 & [cfg.openMM] == 0 & ...
                 [cfg.closeMM] == 0 & ~[cfg.fill] & [cfg.conn] == 8);    % Otsu come in compareModalities

for mo = ["V", "ADC"]
    C = loadStartSlices(T, mo); n = numel(C);

    %% Dice di ogni configurazione su ogni paziente
    D = zeros(numel(cfg), n);
    for g = 1:numel(cfg)
        P = cfg(g); P.polarity = polarity.(char(mo)); P.seedRadiusMM = 3;
        for i = 1:n, D(g, i) = diceCoef(segmentSlice(C(i).img, C(i).seed, C(i).sp, P), C(i).ref); end
        if mod(g, 300) == 0, fprintf('%s: %d/%d configurazioni (%.0f s)\n', mo, g, numel(cfg), toc); end
    end
    med = median(D, 2);
    R = struct2table(cfg); R.dice_median = med; R.dice_mean = mean(D, 2);
    R.dice_q25 = prctile(D, 25, 2); R.dice_q75 = prctile(D, 75, 2);
    writetable(sortrows(R, 'dice_median', 'descend'), fullfile(outDir, "grid_" + mo + ".csv"));

    %% Riepilogo: migliori configurazioni e leave-one-out
    [~, best] = max(med);
    io = find(isOtsu); [~, b] = max(med(io)); bestOtsu = io(b);
    [loo, looCfg] = looSelect(D);
    variants = ["Otsu di partenza"; "Otsu tarato (miglior config Otsu)"; "Miglior config assoluta"; ...
                "Leave-one-out (stima onesta)"];
    X = [D(otsuStart, :); D(bestOtsu, :); D(best, :); loo'];
    cfgTxt = [describeParams(cfg(otsuStart), names); describeParams(cfg(bestOtsu), names); ...
              describeParams(cfg(best), names); sprintf('%d config diverse scelte nei %d fold', numel(unique(looCfg)), n)];
    Sm = table(variants, cfgTxt, median(X, 2), mean(X, 2), prctile(X, 25, 2), prctile(X, 75, 2), ...
        'VariableNames', {'variante', 'configurazione', 'dice_median', 'dice_mean', 'dice_q25', 'dice_q75'});
    writetable(Sm, fullfile(outDir, "summary_" + mo + ".csv"));
    writetable(table([C.id]', X(1, :)', X(2, :)', X(3, :)', loo, 'VariableNames', ...
        {'patient_id', 'dice_otsu_partenza', 'dice_otsu_tarato', 'dice_best', 'dice_loo'}), ...
        fullfile(outDir, "per_patient_" + mo + ".csv"));
    fprintf('\n== %s ==\n', mo); disp(Sm);

    %% Figura: miglior configurazione su ogni paziente
    P = cfg(best); P.polarity = polarity.(char(mo)); P.seedRadiusMM = 3;
    for i = n:-1:1
        [m, box] = segmentSlice(C(i).img, C(i).seed, C(i).sp, P); id = char(C(i).id);
        tiles(i) = struct('img', C(i).img, 'ref', C(i).ref, 'mask', m, 'box', box, 'seeds', C(i).seed, ...
                          'label', sprintf('%s  Dice %.2f', id(end-5:end), D(best, i)));
    end
    saveContactSheet(fullfile(outDir, "contact_best_" + mo + ".png"), tiles, mo + " - " + describeParams(P, names));
    clear tiles
end
fprintf('\nFatto in %.0f s. Risultati in %s\n', toc, outDir);
end
