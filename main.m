% MAIN  Riproduce tutti gli esperimenti in ordine (eseguire dalla cartella gcm_project).
%   Prima vanno preparati i dati con gli script Python (README.md, sezione Data). Tempo totale: circa mezz'ora.
%   Per provare solo il metodo finale: checkPipeline (circa 1 minuto).
addpath('src', 'experiments');
checkPipeline;
checkDataset;                % i casi dev si leggono e sono coerenti
compareModalities;           % quale sequenza RM usare (V, T2_FS, ADC)
tuneStartSlice;              % soglia e morfologia sulla slice di partenza, un seme
tuneMultiSeed;               % un seme per ogni pezzo del tumore
tunePropagation;             % propagazione tra slice: v1 e v2 (metodo finale)
tuneThicknessStop;           % spessore della parete come criterio di arresto (scartato)
ablationRobustness("dev");   % contributo di ogni passo, robustezza, casi di errore
ablationRobustness("test");  % lo stesso sul test, metodo finale senza tarature
analyzeTest;                 % valutazione finale sul test
analyzeSeedFallback;         % metodo finale con e senza seedFallback (risultati vuoti senza la correzione)
renderTestCases3D;           % figura 3D dei casi migliore, mediano e peggiore
plotPipelineSteps;           % figura dei passi del metodo su un paziente dev (relazione)
