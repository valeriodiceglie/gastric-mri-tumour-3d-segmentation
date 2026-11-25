function checkDataset()
%CHECKDATASET Verifica che i casi del set dev si leggano e abbiano dimensioni coerenti.
%   Stampa per ogni caso dimensioni, spaziatura, slice con tumore, slice di area massima e
%   intervallo di intensita'; poi mostra una slice del primo caso con il contorno dell'esperto.
T = readCases('manifest/dev_cases.csv');
fprintf('%d coppie volume+maschera nel set dev\n', height(T));
for i = 1:height(T)
    [I, S, sp] = readCase(T.img(i), T.seg(i));
    area = squeeze(sum(S, [1 2])); [~, k] = max(area);
    fprintf('%s %-5s  size=%s  sp=[%.2f %.2f %.2f]  tumor_slices=%d  slice_max=%d  Imin/max=[%g %g]\n', ...
        T.patient_id(i), T.modality(i), mat2str(size(I)), sp, nnz(area), k, min(I(:)), max(I(:)));
end

%% Controllo visivo: slice di area massima del primo caso
[I, S] = readCase(T.img(1), T.seg(1));
[~, k] = max(squeeze(sum(S, [1 2])));
figure; imshow(mat2gray(I(:, :, k)'), []); hold on;
visboundaries(S(:, :, k)', 'Color', 'r'); title(sprintf('%s %s slice %d', T.patient_id(1), T.modality(1), k));
end
