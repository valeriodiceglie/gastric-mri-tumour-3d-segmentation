% CHECKPIPELINE  Controllo rapido (circa 1 minuto) da eseguire dopo ogni modifica al codice.
%   1) gli helper danno i risultati attesi su esempi piccoli;
%   2) il metodo finale (segmentCase) riproduce il Dice per volume salvato per i 20 pazienti dev
%      (results/ablationRobustness/dev/ablation_per_patient.csv, colonna v6).
%   Eseguire dalla cartella gcm_project:  checkPipeline
addpath('src', 'experiments');

%% 1) Helper
assert(isequal(looSelect([0.9 0.9 0.1; 0.5 0.5 0.8]), [0.5; 0.5; 0.1]));   % il paziente lasciato fuori non sceglie
m = false(4, 4, 3); m(2:3, 2:3, 2) = true;
assert(isequal(volMetrics(m, m), [1 1 0 0]));
s = m; s(:, :, 3) = m(:, :, 2); assert(isequal(volMetrics(m, s), [2/3 0.5 0 0]));
g = paramGrid(struct('a', 0, 'b', ''), 'a', [1 2], 'b', {'x', 'y', 'z'});
assert(isequal([g.a], [1 1 1 2 2 2]) && isequal({g.b}, {'x', 'y', 'z', 'x', 'y', 'z'}));   % l'ultimo campo varia piu' in fretta
B = false(20); B(2:4, 2:4) = true; B(10:17, 10:17) = true; B(19, 19) = true;
assert(isequal(pickSeeds(B, [1 1], 5), [13 13; 3 3]));                    % dal pezzo piu' grande; il pixel isolato e' scartato
assert(isequal(pickSeeds(B, [1 1], 100), [13 13]));                        % il pezzo piu' grande resta sempre
Pc = methodConfig(); Pc.openMM = 0; Pc.sigmaMM = 0;                        % regione chiara con il pixel del seme scuro
I2 = zeros(40); I2(11:30, 11:30) = 1 + 0.2 * mod((11:30)' + (11:30), 2); I2(20, 20) = 0.2;
assert(any(segmentSlice(I2, [20 20], [1 1], Pc), 'all'));                  % con seedFallback la crescita parte
Pc.seedFallback = false; assert(~any(segmentSlice(I2, [20 20], [1 1], Pc), 'all'));   % senza, resta vuota

%% 2) Metodo finale contro i risultati salvati
T = readCases('manifest/dev_cases.csv'); T = T(T.modality == "V", :);
R = readCases('results/ablationRobustness/dev/ablation_per_patient.csv');
for i = 1:height(T)
    [I, S, sp] = readCase(T.img(i), T.seg(i));
    d = diceCoef(segmentCase(int16(I), S, sp), S);
    ref = str2double(R.v6(R.patient_id == T.patient_id(i)));
    assert(abs(d - ref) < 1e-9, '%s: Dice %.6f, salvato %.6f', T.patient_id(i), d, ref);
end
fprintf('checkPipeline: OK (%d pazienti dev identici ai risultati salvati)\n', height(T));
