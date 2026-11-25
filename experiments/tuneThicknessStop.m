function tuneThicknessStop()
%TUNETHICKNESSSTOP Spessore della parete come criterio di arresto della propagazione (set dev, V).
%   Premessa dalla letteratura: il tumore ispessisce la parete gastrica; la parete sana distesa
%   misura 2-5 mm. Quindi una slice segmentata troppo sottile potrebbe non contenere tumore.
%   1) Diagnostica: con il metodo finale, le slice segmentate dove c'e' tumore (TP) sono piu'
%      spesse di quelle segmentate dove non c'e' (FP)?
%   2) Griglia: soglia di spessore x statistica (wallThickness) x arresto relativo, con leave-one-out.
%   Esito: nessun guadagno, il criterio e' scartato.
%   Uscite in results/tuneThicknessStop/: thickness_diagnostic.csv, grid_thickness_V.csv,
%   per_patient_thickness_V.csv
T = readCases('manifest/dev_cases.csv');
outDir = 'results/tuneThicknessStop'; [~, ~] = mkdir(outDir); tic;
[P, Q0] = methodConfig();
Q0.minThickMM = 0; Q0.thickStat = 'max';                       % criterio di spessore spento
C = loadVolumes(T, "V"); n = numel(C);

%% 1) Diagnostica: spessore delle slice segmentate (TP, FP) e della maschera esperta
D = table();
for i = 1:n
    c = C(i); pxmm = mean(c.sp(1:2)); mm2 = prod(c.sp(1:2));
    M3 = propagateVolume(c.vol, c.sp, c.k0, c.seeds, P, Q0);
    for k = 1:size(M3, 3)
        m = M3(:, :, k); e = c.ref(:, :, k);
        if any(m(:))
            if any(e(:)), grp = "slice con tumore (TP)"; else, grp = "slice senza tumore (FP)"; end
            D = [D; thicknessRow(c.id, k, grp, m, pxmm, mm2)]; %#ok<AGROW>
        end
        if any(e(:)), D = [D; thicknessRow(c.id, k, "maschera esperto", e, pxmm, mm2)]; end %#ok<AGROW>
    end
end
writetable(D, fullfile(outDir, 'thickness_diagnostic.csv'));
G = groupsummary(D, 'gruppo', {'median', @(x) prctile(x, 25), @(x) prctile(x, 75)}, ...
                 {'spessore_max_mm', 'spessore_skel_mm', 'area_mm2'});
fprintf('== Diagnostica spessore (metodo finale) ==\n'); disp(G);

%% 2) Griglia del criterio di spessore
grid = paramGrid(Q0, 'minThickMM', [0 4 6 8 10 12], 'thickStat', {'max', 'skel'}, 'minRelArea', [0 0.15]);
grid = grid(~([grid.minThickMM] == 0 & strcmp({grid.thickStat}, 'skel')));   % spessore 0 = spento: una volta sola
nq = numel(grid); DS = zeros(nq, n); EX = DS; RC = DS;
for q = 1:nq
    for i = 1:n
        vm = volMetrics(propagateVolume(C(i).vol, C(i).sp, C(i).k0, C(i).seeds, P, grid(q)), C(i).ref);
        DS(q, i) = vm(1); RC(q, i) = vm(2); EX(q, i) = vm(3);
    end
    fprintf('%d/%d configurazioni (%.0f s)\n', q, nq, toc);
end
R = table(string({grid.thickStat})', [grid.minThickMM]', [grid.minRelArea]', median(DS, 2), mean(DS, 2), ...
    prctile(DS, 25, 2), median(RC, 2), mean(EX, 2), 'VariableNames', {'statistica', 'spessore_min_mm', ...
    'arresto_relativo', 'dice3d_median', 'dice3d_mean', 'dice3d_q25', 'slice_trovate_median', 'slice_in_piu_media'});
writetable(sortrows(R, 'dice3d_median', 'descend'), fullfile(outDir, 'grid_thickness_V.csv'));

%% Confronto con il metodo finale (spessore spento, arresto relativo 0,15)
base = find([grid.minThickMM] == 0 & [grid.minRelArea] == Q0.minRelArea);
[~, best] = max(median(DS, 2)); [loo, looCfg] = looSelect(DS);
fprintf('\n== Griglia spessore (V) ==\n'); disp(sortrows(R, 'dice3d_median', 'descend'));
fprintf('Metodo finale: Dice3D mediana %.3f, slice in piu'' %.2f\n', median(DS(base, :)), mean(EX(base, :)));
fprintf('Migliore: %s %g mm, arresto relativo %.2f -> mediana %.3f, slice in piu'' %.2f\n', grid(best).thickStat, ...
    grid(best).minThickMM, grid(best).minRelArea, median(DS(best, :)), mean(EX(best, :)));
fprintf('Leave-one-out: mediana %.3f (media %.3f), %d config diverse nei %d fold\n', median(loo), mean(loo), ...
    numel(unique(looCfg)), n);
writetable(table([C.id]', DS(base, :)', DS(best, :)', loo, EX(base, :)', EX(best, :)', 'VariableNames', ...
    {'patient_id', 'dice3d_metodo_finale', 'dice3d_best', 'dice3d_loo', 'slice_in_piu_finale', 'slice_in_piu_best'}), ...
    fullfile(outDir, 'per_patient_thickness_V.csv'));
fprintf('\nFatto in %.0f s.\n', toc);
end

function r = thicknessRow(id, k, grp, m, pxmm, mm2)
    r = table(id, k, grp, wallThickness(m, pxmm, 'max'), wallThickness(m, pxmm, 'skel'), nnz(m) * mm2, ...
        'VariableNames', {'patient_id', 'slice', 'gruppo', 'spessore_max_mm', 'spessore_skel_mm', 'area_mm2'});
end
