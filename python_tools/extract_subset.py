#!/usr/bin/env python3
"""
extract_subset.py  Estrazione dallo zip di un sottoinsieme di pazienti e modalita'.
Copia dallo zip SOLO i file richiesti, con nomi normalizzati:
   data/<split>/<paziente>/<MODALITA>_img.nii.gz   e   <MODALITA>_seg.nii.gz
Scrive anche manifest/<split>_cases.csv con i percorsi relativi (letto dagli script MATLAB).
Uso: python extract_subset.py <GCM2025.zip> <cartella_progetto> <split> <modalita,separate,da,virgola> [file_split colonna da a]
Esempio: python extract_subset.py ../../GCM2025.zip .. dev V,T2_FS,ADC
Set test (474 pazienti), a blocchi di pazienti:
         python extract_subset.py ../GCM2025.zip . test V manifest/split_final.csv split_final 0 120
"""
import sys, os, zipfile, pandas as pd
ZIP, PROJ, SPLIT, MODS = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4].split(',')
m = pd.read_csv(f'{PROJ}/manifest/manifest_modalities.csv', dtype={'patient_id': str})
SFILE = sys.argv[5] if len(sys.argv) > 5 else f'{PROJ}/manifest/split.csv'
SCOL = sys.argv[6] if len(sys.argv) > 6 else 'split'
s = pd.read_csv(SFILE if SFILE.startswith('/') or SFILE.startswith('.') else f'{PROJ}/{SFILE}', dtype={'patient_id': str})
ids = sorted(s[s[SCOL] == SPLIT].patient_id.tolist())
ALL_IDS = ids
if len(sys.argv) > 8: ids = ids[int(sys.argv[7]):int(sys.argv[8])]
sel = m[m.patient_id.isin(ids) & m.modality.isin(MODS)]
assert (sel.status == 'ok').all(), 'ci sono modalita\' non regolari nel sottoinsieme'
zf = zipfile.ZipFile(ZIP); rows = []; tot = 0
for r in sel.itertuples():
    d = f'{PROJ}/data/{SPLIT}/{r.patient_id}'; os.makedirs(d, exist_ok=True)
    for kind, member in (('img', r.zip_img), ('seg', r.zip_seg)):
        dst = f'{d}/{r.modality}_{kind}.nii.gz'
        if not os.path.exists(dst) or os.path.getsize(dst) != zf.getinfo(member).file_size:
            data = zf.read(member); open(dst, 'wb').write(data)
        tot += os.path.getsize(dst)
    rows.append(dict(patient_id=r.patient_id, modality=r.modality,
                     img=f'data/{SPLIT}/{r.patient_id}/{r.modality}_img.nii.gz',
                     seg=f'data/{SPLIT}/{r.patient_id}/{r.modality}_seg.nii.gz'))
# il file dei casi elenca sempre TUTTI i pazienti dello split (anche se si estrae a blocchi)
c = pd.DataFrame([dict(patient_id=p, modality=mo, img=f'data/{SPLIT}/{p}/{mo}_img.nii.gz', seg=f'data/{SPLIT}/{p}/{mo}_seg.nii.gz')
                  for p in ALL_IDS for mo in MODS]).sort_values(['patient_id', 'modality'])
c.to_csv(f'{PROJ}/manifest/{SPLIT}_cases.csv', index=False)
print(f'{len(ids)} pazienti, {len(rows)} coppie volume+maschera, {tot/1e6:.0f} MB in data/{SPLIT}/')
