function v = volMetrics(M3, S3)
%VOLMETRICS Metriche per volume: [Dice 3D, frazione di slice del tumore trovate,
%   slice segmentate senza tumore, arrivo alla prima o all'ultima slice (0/1)].
    pe = squeeze(any(M3, [1 2])); ex = squeeze(any(S3, [1 2]));
    v = [diceCoef(M3, S3), mean(pe(ex)), nnz(pe & ~ex), double(pe(1) || pe(end))];
end
