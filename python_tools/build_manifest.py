#!/usr/bin/env python3
"""
build_manifest.py  Manifest dei 500 pazienti (primo passo della preparazione dei dati).
Legge GCM2025.zip SENZA estrarlo e produce due CSV in manifest/:
  - manifest_modalities.csv : una riga per (paziente, modalita') con stato, anomalie,
                              header (dimensioni, spaziatura, orientamento) e statistiche della maschera
  - manifest_patients.csv   : una riga per paziente, con i dati clinici dell'Excel e
                              il flag di idoneita' (V, T2_FS e ADC tutte regolari)
Uso: python build_manifest.py <GCM2025.zip> <cartella_manifest>
"""
import sys, gzip, zipfile, collections, re
import numpy as np, pandas as pd, nibabel as nib

ZIP, OUT = sys.argv[1], sys.argv[2]
CORE = ['V', 'T2_FS', 'ADC']          # modalita' presenti in tutti i 500 pazienti
zf = zipfile.ZipFile(ZIP)

# ---- 1. raggruppa i file per paziente/modalita' ------------------------------------------
entries = collections.defaultdict(lambda: collections.defaultdict(list))
nested = set()
for zi in zf.infolist():
    if zi.is_dir():
        continue
    parts = zi.filename.split('/')
    name = parts[-1]
    if parts[0] != 'data' or '__MACOSX' in parts or name.startswith('._'):
        continue                        # file di sistema macOS: ignorati
    if len(parts) == 4:                 # data/<paziente>/<modalita'>/<file>
        entries[parts[1]][parts[2]].append(zi)
    elif len(parts) == 5:               # data/<A>/<paziente>/<modalita'>/<file>: copia annidata
        nested.add(parts[2])

def resolve(pid, zis):
    """Stato della coppia (volume, maschera): ok | ok_extra | irregular."""
    names = {z.filename.split('/')[-1]: z for z in zis}
    img_n, seg_n = f'{pid}.nii.gz', f'{pid}seg.nii.gz'
    if img_n in names and seg_n in names:
        extras = sorted(n for n in names if n not in (img_n, seg_n))
        if not extras:
            return 'ok', names[img_n], names[seg_n], ''
        return 'ok_extra', names[img_n], names[seg_n], 'file_extra:' + '|'.join(extras)
    return 'irregular', None, None, 'trovati:' + '|'.join(sorted(names))

def read_header(zi):
    with zf.open(zi) as f, gzip.GzipFile(fileobj=f) as g:
        return nib.Nifti1Header(binaryblock=g.read(348))

def read_seg(zi):
    raw = gzip.decompress(zf.read(zi))
    img = nib.Nifti1Image.from_bytes(raw)
    return np.asanyarray(img.dataobj)

rows = []
for n, pid in enumerate(sorted(entries)):
    for mod, zis in sorted(entries[pid].items()):
        status, zi_img, zi_seg, issue = resolve(pid, zis)
        r = dict(patient_id=pid, modality=mod, status=status, issues=issue,
                 zip_img=zi_img.filename if zi_img else '', zip_seg=zi_seg.filename if zi_seg else '')
        if pid in nested:
            r['issues'] = (r['issues'] + ';' if r['issues'] else '') + 'esiste_copia_annidata'
        if zi_img and mod in CORE:
            try:
                h, hs = read_header(zi_img), read_header(zi_seg)
                shp, sp = h.get_data_shape(), h.get_zooms()
                seg = read_seg(zi_seg)
                m = seg > 0
                area = m.sum(axis=(0, 1))
                ks = np.nonzero(area)[0]
                r.update(shape_x=shp[0], shape_y=shp[1], shape_z=shp[2] if len(shp) > 2 else 1,
                         sp_x=round(float(sp[0]), 4), sp_y=round(float(sp[1]), 4), sp_z=round(float(sp[2]), 4),
                         orient=''.join(nib.aff2axcodes(h.get_best_affine())),
                         same_grid=bool(tuple(hs.get_data_shape()[:3]) == tuple(shp[:3])),
                         seg_labels='|'.join(str(int(v)) for v in np.unique(seg)),
                         tumor_voxels=int(m.sum()), tumor_slices=int(len(ks)))
                if len(ks):
                    kmax = int(np.argmax(area))
                    c0, c1 = np.nonzero(m[:, :, kmax])
                    a_mm2 = float(area[kmax]) * sp[0] * sp[1]
                    r.update(kmax=kmax, max_area_px=int(area[kmax]), centroid_x=round(float(c0.mean()), 2),
                             centroid_y=round(float(c1.mean()), 2),
                             eq_diam_mm=round(float(np.sqrt(4 * a_mm2 / np.pi)), 2),
                             volume_mm3=round(float(m.sum()) * sp[0] * sp[1] * sp[2], 1))
            except Exception as e:                      # file illeggibile: lo si segnala
                r['status'] = 'irregular'
                r['issues'] = (r['issues'] + ';' if r['issues'] else '') + f'errore_lettura:{type(e).__name__}'
        rows.append(r)
    if n % 50 == 0:
        print(f'{n}/{len(entries)} pazienti', file=sys.stderr, flush=True)

mod = pd.DataFrame(rows)

# ---- 2. idoneita' per paziente -----------------------------------------------------------
def eligible(g):
    d = g[g.modality.isin(CORE)].set_index('modality')
    if not all(m in d.index for m in CORE):
        return False
    for m in CORE:
        x = d.loc[m]
        if x.status != 'ok' or not bool(x.get('same_grid', False)):
            return False
        if str(x.get('seg_labels', '')) not in ('0|1', '1'):      # maschera binaria
            return False
        if not (x.get('tumor_voxels', 0) > 0):
            return False
    return True

pat = pd.DataFrame({'patient_id': sorted(entries)})
pat['eligible'] = [eligible(mod[mod.patient_id == p]) for p in pat.patient_id]
pat['n_modalities'] = [mod[mod.patient_id == p].shape[0] for p in pat.patient_id]
pat['has_nested_copy'] = pat.patient_id.isin(nested)

# ---- 3. dati clinici dall'Excel (intestazioni in cinese, rinominate per posizione) --------
with zf.open('GCM2025.xlsx') as f:
    x = pd.read_excel(f)
x.columns = ['patient_id_int', 'mri_date', 'CEA_ng_ml', 'CA125_U_ml', 'CA19_9_U_ml', 'max_thickness_mm']
def num(s):
    s = pd.to_numeric(s, errors='coerce')
    return s
x['patient_id'] = x.patient_id_int.astype(int).map(lambda v: f'{v:010d}')
for c in ['CEA_ng_ml', 'CA125_U_ml', 'CA19_9_U_ml', 'max_thickness_mm']:
    x[c] = num(x[c])
pat = pat.merge(x.drop(columns='patient_id_int'), on='patient_id', how='left')

mod.to_csv(f'{OUT}/manifest_modalities.csv', index=False)
pat.to_csv(f'{OUT}/manifest_patients.csv', index=False)

# ---- 4. riepilogo a video -------------------------------------------------------------------
print('\n== RIEPILOGO ==')
print('pazienti:', len(pat), '| idonei (V, T2_FS e ADC regolari):', int(pat.eligible.sum()))
print('stato per modalita\':'); print(mod.groupby(['modality', 'status']).size().unstack(fill_value=0).to_string())
print('\nanomalie (stato diverso da ok):')
print(mod[mod.status != 'ok'][['patient_id', 'modality', 'status', 'issues']].to_string(index=False))
print('\nmatching Excel:', int(pat.mri_date.notna().sum()), 'pazienti su', len(pat))
