#!/usr/bin/env python3
"""
make_split.py  Split iniziale dev/test stratificato per volume del tumore.
Assegna i pazienti idonei a: dev (20), test (40), unused (gli altri).
- stratificato per terzili del volume del tumore (maschera V), per avere tumori piccoli, medi e grandi
  in ogni gruppo;
- seed fisso: lo split e' riproducibile;
- il gruppo test non si guarda durante la taratura: serve solo per la valutazione finale.
  (Lo split usato alla fine, con test = tutti gli altri idonei, e' manifest/split_final.csv.)
Uso: python make_split.py <cartella_manifest> [seed]
"""
import sys, numpy as np, pandas as pd
OUT = sys.argv[1]; SEED = int(sys.argv[2]) if len(sys.argv) > 2 else 42
N_DEV, N_TEST = 20, 40
m = pd.read_csv(f'{OUT}/manifest_modalities.csv', dtype={'patient_id': str})
p = pd.read_csv(f'{OUT}/manifest_patients.csv', dtype={'patient_id': str})
vol = m[m.modality == 'V'].set_index('patient_id').volume_mm3
p['volume_V_mm3'] = p.patient_id.map(vol)
el = p[p.eligible].copy()
el['tercile'] = pd.qcut(el.volume_V_mm3, 3, labels=['piccolo', 'medio', 'grande']).astype(str)
rng = np.random.default_rng(SEED)
split = {pid: 'unused' for pid in p.patient_id}
def alloc(n):                       # ripartisce n il piu' equamente possibile sui 3 terzili
    base = [n // 3] * 3
    for i in range(n - sum(base)): base[i] += 1
    return dict(zip(['piccolo', 'medio', 'grande'], base))
dev_n, test_n = alloc(N_DEV), alloc(N_TEST)
for t in ['piccolo', 'medio', 'grande']:
    ids = el[el.tercile == t].patient_id.to_numpy()
    pick = rng.permutation(ids)
    for pid in pick[:dev_n[t]]: split[pid] = 'dev'
    for pid in pick[dev_n[t]:dev_n[t] + test_n[t]]: split[pid] = 'test'
p['split'] = p.patient_id.map(split)
p = p.merge(el[['patient_id', 'tercile']], on='patient_id', how='left')
p[['patient_id', 'eligible', 'split', 'tercile', 'volume_V_mm3', 'max_thickness_mm']].to_csv(f'{OUT}/split.csv', index=False)
print(p.groupby(['split', 'tercile'], dropna=False).size().unstack(fill_value=0).to_string())
print('\nvolume V (cm3) mediano per gruppo:')
print((p.groupby('split').volume_V_mm3.median() / 1000).round(1).to_string())
