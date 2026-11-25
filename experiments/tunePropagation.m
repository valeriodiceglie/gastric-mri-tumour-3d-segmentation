function tunePropagation()
%TUNEPROPAGATION Taratura della propagazione tra slice e Dice per volume (set dev, V e ADC).
%   La slice di partenza (area massima), i semi (uno per pezzo per V, uno solo per ADC) e la
%   segmentazione di una slice sono fissi (methodConfig); si tarano solo i parametri della
%   propagazione, con griglia e leave-one-out. Due versioni:
%     v1  semi o maschera precedente, vincolo di distanza, arresto su area e sovrapposizione.
%         Volumi ritagliati alle slice del tumore +- 8: le fughe oltre il ritaglio non sono
%         contate, quindi il Dice e' ottimistico (errore trovato dopo, corretto in v2).
%     v2  in piu' limite di crescita, ROI che segue la maschera, arresto relativo all'area
%         massima; volume intero. Da qui viene la configurazione finale (methodConfig).
%   Due riferimenti: solo la slice di partenza; semi dall'esperto su ogni slice (stessa
%   segmentazione, ma ogni slice riceve i semi giusti: limite superiore della propagazione).
%   Uscite in results/tunePropagation/<v1|v2>/: grid_, summary_, per_patient_<mod>.csv e
%   area_profiles_<mod>.png (area per slice, esperto contro propagazione).
T = readCases('manifest/dev_cases.csv'); tic;
[~, Q, Qv1] = methodConfig();
vary.v1 = {'marker', {'points', 'mask'}, 'erodeMM', [0 1 2], 'clipMM', [0 2 4 8], 'minAreaMM2', [20 50], ...
           'minOverlap', [0 0.2 0.4]};
vary.v2 = {'marker', {'points', 'mask'}, 'clipMM', [1 2 4], 'maxGrowth', [1.2 1.5 Inf], 'roiFollow', [false true], ...
           'minRelArea', [0 0.15 0.3]};
base.v1 = Qv1; base.v2 = Q;
zMargin.v1 = 8; zMargin.v2 = Inf;
maxSeeds.V = Inf; maxSeeds.ADC = 1;
for v = ["v1", "v2"]
    outDir = fullfile('results', 'tunePropagation', v); [~, ~] = mkdir(outDir);
    for mo = ["V", "ADC"]
        b = base.(v); b.maxSeeds = maxSeeds.(mo);
        tuneOne(T, mo, paramGrid(b, vary.(v){:}), vary.(v)(1:2:end), zMargin.(v), outDir);
        fprintf('%s %s fatto (%.0f s)\n', v, mo, toc);
    end
end
end

function tuneOne(T, mo, grid, names, zMargin, outDir)
% Griglia della propagazione per una modalita': CSV e figura in outDir.
    P = methodConfig(mo);
    C = loadVolumes(T, mo, zMargin); n = numel(C);
    for i = 1:n, C(i).seeds = C(i).seeds(1:min(end, grid(1).maxSeeds), :); end

    %% Riferimenti: solo slice di partenza; semi dall'esperto su ogni slice
    mStart = zeros(n, 4); mRef = zeros(n, 4);
    for i = 1:n
        c = C(i);
        M3 = false(size(c.ref));
        M3(:, :, c.k0) = segmentSlice(double(c.vol(:, :, c.k0)), c.seeds, c.sp(1:2), P);
        mStart(i, :) = volMetrics(M3, c.ref);
        M3 = false(size(c.ref));
        for k = find(squeeze(any(c.ref, [1 2])))'
            sd = pickSeeds(c.ref(:, :, k), c.sp(1:2), grid(1).minSeedAreaMM2);
            M3(:, :, k) = segmentSlice(double(c.vol(:, :, k)), sd(1:min(end, grid(1).maxSeeds), :), c.sp(1:2), P);
        end
        mRef(i, :) = volMetrics(M3, c.ref);
    end

    %% Griglia: Dice per volume, slice del tumore trovate, slice in piu', arrivo al bordo
    nq = numel(grid); D = zeros(nq, n); REC = D; EXTRA = D; HIT = D;
    for q = 1:nq
        for i = 1:n
            vm = volMetrics(propagateVolume(C(i).vol, C(i).sp, C(i).k0, C(i).seeds, P, grid(q)), C(i).ref);
            D(q, i) = vm(1); REC(q, i) = vm(2); EXTRA(q, i) = vm(3); HIT(q, i) = vm(4);
        end
    end
    med = median(D, 2); [~, best] = max(med);
    [loo, looCfg] = looSelect(D);
    R = struct2table(grid); R.dice3d_median = med; R.dice3d_mean = mean(D, 2);
    R.recall_median = median(REC, 2); R.extra_mean = mean(EXTRA, 2); R.hit_bordo = sum(HIT, 2);
    writetable(sortrows(R, 'dice3d_median', 'descend'), fullfile(outDir, "grid_" + mo + ".csv"));

    %% Riepilogo e risultati per paziente
    variants = ["Solo slice di partenza"; "Propagazione, config migliore"; "Leave-one-out della propagazione"; ...
                "Riferimento: semi dall'esperto su ogni slice"];
    X = [mStart(:, 1)'; D(best, :); loo'; mRef(:, 1)'];
    recall = [median(mStart(:, 2)); median(REC(best, :)); NaN; median(mRef(:, 2))];
    extra = [mean(mStart(:, 3)); mean(EXTRA(best, :)); NaN; mean(mRef(:, 3))];
    cfgTxt = ["-"; describeParams(grid(best), names); sprintf('%d config diverse nei %d fold', numel(unique(looCfg)), n); "-"];
    Sm = table(variants, cfgTxt, median(X, 2), mean(X, 2), prctile(X, 25, 2), prctile(X, 75, 2), recall, extra, ...
        'VariableNames', {'variante', 'configurazione', 'dice3d_median', 'dice3d_mean', 'dice3d_q25', 'dice3d_q75', ...
                          'slice_tumore_trovate_mediana', 'slice_in_piu_media'});
    writetable(Sm, fullfile(outDir, "summary_" + mo + ".csv"));
    writetable(table([C.id]', [C.nTumorSlices]', mStart(:, 1), D(best, :)', loo, mRef(:, 1), REC(best, :)', ...
        EXTRA(best, :)', HIT(best, :)', 'VariableNames', {'patient_id', 'slice_tumore', 'dice3d_solo_partenza', ...
        'dice3d_best', 'dice3d_loo', 'dice3d_riferimento', 'slice_trovate_best', 'slice_in_piu_best', 'hit_bordo_best'}), ...
        fullfile(outDir, "per_patient_" + mo + ".csv"));
    fprintf('\n== %s ==\n', mo); disp(Sm);

    %% Figura: area per slice, esperto (rosso) e propagazione con la config migliore (verde)
    f = figure('Visible', 'off', 'Position', [0 0 1400 1000]); theme(f, 'light');
    tl = tiledlayout(f, 4, 5, 'TileSpacing', 'compact');
    for i = 1:n
        c = C(i); mm2 = c.sp(1) * c.sp(2);
        M3 = propagateVolume(c.vol, c.sp, c.k0, c.seeds, P, grid(best));
        ae = squeeze(sum(c.ref, [1 2])) * mm2; ap = squeeze(sum(M3, [1 2])) * mm2;
        nexttile(tl); plot(ae, 'r-', 'LineWidth', 1.2); hold on; plot(ap, 'g-', 'LineWidth', 1.2); xline(c.k0, 'k:');
        ks = find(ae > 0 | ap > 0); xlim([max(min(ks) - 3, 1), min(max(ks) + 3, numel(ae))]);
        id = char(c.id); title(sprintf('%s  Dice3D %.2f', id(end-5:end), D(best, i)), 'FontSize', 8);
        if i == 1, ylabel('area (mm^2)'); legend({'expert', 'propagation', 'start'}, 'FontSize', 6, 'Location', 'best'); end
    end
    title(tl, sprintf('%s - area per slice (%s)', mo, describeParams(grid(best), names)), 'FontSize', 9, 'Interpreter', 'none');
    exportgraphics(f, fullfile(outDir, "area_profiles_" + mo + ".png"), 'Resolution', 100); close(f);
end
