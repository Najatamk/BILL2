#!/usr/bin/env python3
"""
Unified KHV Pangenome & Metagenomic Analysis Pipeline
====================================================
Intègre toutes les fonctionnalités de `main.py` et `analysis_pipeline.py` dans un
pipeline modulaire exécutable en une seule commande ou par étapes.

Principales capacités:
  1. Pré-QC (statistiques assembly)
  2. Kraken2 classification + décontamination (wrapper + SLURM option)
  3. Post-QC + comparaison pré/post
  4. QUAST (pré et post si FASTA nettoyés)
  5. BUSCO (pré et post)
  6. Construction pangenome (stubs exécutables): Minigraph-Cactus & PGGB
  7. Variants (normalisation, filtration, comparaison)
  8. Analyses avancées: phylogénie, syntenie, impact fonctionnel, sélection
  9. Benchmarking & monitoring
 10. Résumés globaux + visualisations
 11. Système de validation & reproductibilité

Usage rapide:
  python unified_pipeline.py --genomes p15:data/p15_khv.fasta p90:data/p90_khv.fasta --all
  python unified_pipeline.py --genomes p15:... p90:... --steps pre-qc,kraken,post-qc,quast,busco,summary

Chaque étape est idempotente avec détection de cache quand possible.
"""
from __future__ import annotations
import os, sys, argparse, subprocess, shutil, json, csv, gzip, statistics, re, textwrap, time
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Callable
from dataclasses import dataclass, field

# Imports optionnels pour plots
try:
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt
except ImportError:
    pd = None
    sns = None
    plt = None

# ===================== Configuration Portable =====================
# Détection automatique des chemins d'outils via variable d'environnement
# ou recherche dans des emplacements standards
TOOLS_ROOT = os.environ.get('PIPELINE_TOOLS_DIR', '/data/tools')

# AJOUT: répertoires supplémentaires de recherche d'outils
# Ces chemins sont testés s'ils existent, pas d'erreur si absents
EXTRA_SEARCH_DIRS = [
    Path(TOOLS_ROOT) / 'quast',
    Path(TOOLS_ROOT) / 'kraken2',
    Path(TOOLS_ROOT) / 'pggb',
    Path(TOOLS_ROOT) / 'cactus',
    Path(TOOLS_ROOT) / 'snakemake',
    Path(TOOLS_ROOT) / 'bcftools-1.21',
    Path(TOOLS_ROOT) / 'bowtie2',
    Path(TOOLS_ROOT) / 'metaeuk',  # Ajout pour la dépendance de BUSCO
]

# ===================== Configuration & Constantes =====================
RESULTS = Path('results')
PLOTS = RESULTS / 'plots'
REPORTS = RESULTS / 'reports'
KRAKEN_DIR = RESULTS / 'kraken'
PANGENOME_DIR = RESULTS / 'pangenome'
VALIDATION_DIR = RESULTS / 'validation'

# Base de données Kraken2 - configurable via variable d'environnement
DB_DEFAULT = os.environ.get('KRAKEN2_DB', str(Path(TOOLS_ROOT) / 'kraken2' / 'kraken2_db'))

EXPECTED_DEFAULT_SIZE = 295146
EXPECTED_DEFAULT_GC = 59.2

# ===================== Utilitaires Log =====================
COLORS = {
    'reset':'\033[0m','blue':'\033[94m','green':'\033[92m','yellow':'\033[93m','red':'\033[91m','cyan':'\033[96m'
}
def log(msg: str):
    print(f"{COLORS['green']}[PIPELINE]{COLORS['reset']} {msg}")
def warn(msg: str):
    print(f"{COLORS['yellow']}[WARN]{COLORS['reset']} {msg}")
def error(msg: str):
    print(f"{COLORS['red']}[ERROR]{COLORS['reset']} {msg}")

def run_cmd(cmd: List[str], cwd: Optional[Path]=None, capture=False) -> subprocess.CompletedProcess:
    log('CMD: ' + ' '.join(cmd))
    return subprocess.run(cmd, cwd=cwd, text=True, capture_output=capture)

def enrich_path():
    """
    Ajoute automatiquement les répertoires d'outils au PATH.
    Utilise PIPELINE_TOOLS_DIR si défini, sinon essaie des emplacements standards.
    """
    root = Path(TOOLS_ROOT)
    to_add = []
    
    # Dossiers explicitement listés
    for d in EXTRA_SEARCH_DIRS:
        if d.exists():
            to_add.append(d)
            b = d / 'bin'
            if b.exists():
                to_add.append(b)
    
    # Parcours du répertoire racine des outils (profondeur 1)
    if root.exists():
        for sub in root.iterdir():
            if sub.is_dir():
                to_add.append(sub)
                b = sub / 'bin'
                if b.exists():
                    to_add.append(b)
    
    current = os.environ.get('PATH', '')
    added = []
    for p in dict.fromkeys(str(x) for x in to_add):  # ordre + unique
        if p and p not in current:
            current = p + ':' + current
            added.append(p)
    os.environ['PATH'] = current
    if added:
        log(f'Added to PATH: {", ".join(added)}')

# ===================== Data Models =====================
@dataclass
class PipelineConfig:
    genomes: Dict[str, Path]
    reference_key: Optional[str] = None
    pangenome_reference_key: Optional[str] = None
    pangenome_reference_fasta: Optional[str] = None  # Path to external reference FASTA
    kraken_db: str = DB_DEFAULT
    target_taxids: List[int] = field(default_factory=list)
    keep_unclassified: bool = False
    busco_lineage: str = 'viruses_odb10'
    busco_auto_lineage: bool = False  # <--- AJOUT
    quast_threads: int = int(os.environ.get('SLURM_CPUS_PER_TASK', '4'))
    busco_threads: int = int(os.environ.get('SLURM_CPUS_PER_TASK', '4'))
    expected_size: int = EXPECTED_DEFAULT_SIZE
    expected_gc: float = EXPECTED_DEFAULT_GC
    size_tolerance_pct: float = 12.0
    gc_tolerance: float = 5.0
    force_quast: bool = False
    force_busco: bool = False
    skip_kraken_if_present: bool = True
    plots: bool = True
    verbose: bool = False
    steps: List[str] = field(default_factory=list)
    # New: maximum seconds to wait for cleaned FASTA (async Kraken scenario)
    wait_clean_secs: int = 0
    # Pangenome params
    pangenome_threads: int = int(os.environ.get('SLURM_CPUS_PER_TASK', '8'))
    cactus_pipeline_dir: str = 'cactus-snakemake'
    cactus_output_prefix: str = 'khv'
    cactus_config_path: Optional[str] = None
    cactus_skip_if_exists: bool = True
    pggb_segment: int = 5000
    pggb_percent_id: int = 70
    pggb_passes: int = 10
    pggb_threads: int = int(os.environ.get('SLURM_CPUS_PER_TASK', '8'))
    pggb_extra: str = ''
    pggb_skip_if_exists: bool = True

# ===================== Helpers =====================

def ensure_dirs():
    for d in [RESULTS, PLOTS, REPORTS, KRAKEN_DIR, PANGENOME_DIR, VALIDATION_DIR]:
        d.mkdir(parents=True, exist_ok=True)

# FASTA iterator (supports gz)

def fasta_iter(path: Path):
    opener = gzip.open if str(path).endswith(('.gz','.gzip')) else open
    name=None; seq=[]
    with opener(path, 'rt') as fh:
        for line in fh:
            if not line.strip():
                continue
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(seq)
                name = line[1:].strip().split()[0]
                seq=[]
            else:
                seq.append(line.strip())
        if name is not None:
            yield name, ''.join(seq)

def n50(lengths: List[int]) -> int:
    if not lengths: return 0
    s=sorted(lengths, reverse=True)
    half=sum(lengths)/2; acc=0
    for L in s:
        acc+=L
        if acc>=half: return L
    return s[-1]

# ===================== Dépendances =====================
REQUIRED_TOOL_GROUPS = {
    'basic': ['python','bash'],
    'quality': ['quast.py','busco'],
    'contamination': ['kraken2'],
    'variants': ['bcftools','samtools'],
    'pangenome_cactus': ['cactus','singularity'],
    'pangenome_pggb': ['pggb','odgi','vg'],
    'phylogeny': ['fasttree','iqtree'],
}

def which(tool: str) -> bool:
    return shutil.which(tool) is not None

def find_quast():
    """Retourne le binaire quast (quast.py ou quast) si trouvé."""
    cand = shutil.which('quast.py') or shutil.which('quast')
    if cand:
        return cand
    # recherche manuelle
    for base in EXTRA_SEARCH_DIRS:
        q1 = base / 'quast.py'
        q2 = base / 'quast' / 'quast.py'
        if q1.exists():
            return str(q1)
        if q2.exists():
            return str(q2)
    return None

def check_dependencies() -> Dict[str, Dict[str,bool]]:
    log('Checking system dependencies...')
    status={}
    for group, tools in REQUIRED_TOOL_GROUPS.items():
        status[group]={t: which(t) for t in tools}
        missing=[t for t,a in status[group].items() if not a]
        if missing:
            warn(f"Group {group}: missing {', '.join(missing)}")
        else:
            log(f"Group {group}: OK")
    return status

# ===================== QC Pré =====================

def run_pre_qc(cfg: PipelineConfig) -> Path:
    ensure_dirs()
    rows=[]
    for sample, path in cfg.genomes.items():
        if not path.exists():
            warn(f"Missing genome: {sample} -> {path}")
            continue
        lengths=[]; gcs=[]; total_bp=0
        for cid, seq in fasta_iter(path):
            L=len(seq); lengths.append(L); total_bp+=L
            gc=sum(1 for b in seq.upper() if b in 'GC')
            atgc=sum(1 for b in seq.upper() if b in 'ATGC')
            gcs.append((gc/atgc*100) if atgc else 0)
        rows.append({
            'sample': sample,
            'n_contigs': len(lengths),
            'total_bp': total_bp,
            'N50': n50(lengths),
            'mean_len': statistics.mean(lengths) if lengths else 0,
            'median_len': statistics.median(lengths) if lengths else 0,
            'GC_mean': statistics.mean(gcs) if gcs else 0.0
        })
    out = REPORTS / 'pre_decontam_qc_metrics.tsv'
    header = rows[0].keys() if rows else ['sample','n_contigs','total_bp','N50','mean_len','median_len','GC_mean']
    with out.open('w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=header, delimiter='\t')
        writer.writeheader()
        for r in rows:
            writer.writerow(r)
    log(f"Pre-QC metrics -> {out}")
    return out

# ===================== Kraken2 (soumission & wrapper) =====================

def kraken_db_valid(db_path: str) -> bool:
    if not db_path:
        return False
    required = ['hash.k2d', 'opts.k2d', 'taxo.k2d']
    p = Path(db_path)
    return all((p / r).exists() for r in required)

def kraken_outputs_present() -> bool:
    return (KRAKEN_DIR / 'clean' / 'summary_samples.tsv').exists()

def submit_kraken_jobs(cfg: PipelineConfig) -> Optional[Path]:
    if not kraken_db_valid(cfg.kraken_db):
        warn(f"Kraken DB invalide ou incomplète: {cfg.kraken_db} (hash.k2d/opts.k2d/taxo.k2d requis) -> étape sautée")
        return None
    if cfg.skip_kraken_if_present and kraken_outputs_present():
        log('Kraken outputs detected -> skip step')
        return KRAKEN_DIR / 'clean' / 'summary_samples.tsv'
    KRAKEN_DIR.mkdir(parents=True, exist_ok=True)
    # Simple sequential execution (placeholder). For cluster, user can adapt sbatch command
    summary_rows=[]
    for sample, path in cfg.genomes.items():
        if not path.exists():
            warn(f"Skip Kraken (missing): {sample}")
            continue
        # Command (classification only - extraction logic external or stub)
        out_file = KRAKEN_DIR / f"{sample}_kraken.out"
        rep_file = KRAKEN_DIR / f"{sample}_kraken_report.txt"
        cmd=["kraken2","--db", cfg.kraken_db,"--threads","4","--output", str(out_file),"--report", str(rep_file), str(path)]
        r=run_cmd(cmd, capture=True)
        if r.returncode!=0:
            warn(f"Kraken failed for {sample}: {r.stderr.splitlines()[:2] if r.stderr else ''}")
            continue
        # Parse classification output to build mapping contig -> taxid
        classified=0; unclassified=0
        contig_taxid: Dict[str, Optional[int]] = {}
        if out_file.exists():
            for line in out_file.read_text().splitlines():
                if not line: continue
                parts=line.split('\t')
                if len(parts) < 3: continue
                status, read_id, taxid_s = parts[0], parts[1], parts[2]
                if status=='C':
                    classified+=1
                    try:
                        contig_taxid[read_id]=int(taxid_s)
                    except ValueError:
                        contig_taxid[read_id]=None
                elif status=='U':
                    unclassified+=1
                    contig_taxid[read_id]=None
        total_reads = classified + unclassified

        # Determine which taxids to retain
        retain_taxids: List[int] = []
        if cfg.target_taxids:
            retain_taxids = cfg.target_taxids
        else:
            # pick most abundant taxid if none provided
            freq: Dict[int,int] = {}
            for rid, tx in contig_taxid.items():
                if tx is None: continue
                freq[tx]=freq.get(tx,0)+1
            if freq:
                # choose top taxid
                retain_taxids=[sorted(freq.items(), key=lambda x:x[1], reverse=True)[0][0]]
        # Load original FASTA sequences
        seqs: Dict[str,str] = {cid:seq for cid,seq in fasta_iter(path)}
        kept_ids=set()
        for rid, tx in contig_taxid.items():
            if tx is not None and retain_taxids and tx in retain_taxids:
                kept_ids.add(rid)
            elif tx is None and cfg.keep_unclassified:
                kept_ids.add(rid)
        # If filtering eliminates everything, fall back to all sequences
        if not kept_ids:
            warn(f"Filtering for {sample} retained 0 contigs -> keeping all (fallback)")
            kept_ids=set(seqs.keys())
        # Write cleaned FASTA
        clean_dir = KRAKEN_DIR / 'clean'
        clean_dir.mkdir(exist_ok=True)
        clean_fa = clean_dir / f"{sample}_clean.fasta"
        with clean_fa.open('w') as fhc:
            for cid in kept_ids:
                if cid not in seqs: # contig id mismatch safety
                    continue
                fhc.write(f">{cid}\n{seqs[cid]}\n")
        removed = len(seqs) - len(kept_ids)
        kept_percent = (len(kept_ids)/len(seqs)*100) if seqs else 0
        summary_rows.append({
            'sample': sample,
            'total_contigs': len(seqs),
            'classified': classified,
            'unclassified': unclassified,
            'retained_contigs': len(kept_ids),
            'removed_contigs': removed,
            'retained_percent': f"{kept_percent:.2f}",
            'retained_taxids': ';'.join(map(str,retain_taxids)) if retain_taxids else 'NA'
        })
    # Write summary
    if summary_rows:
        clean_dir = KRAKEN_DIR / 'clean'
        clean_dir.mkdir(exist_ok=True)
        summ = clean_dir / 'summary_samples.tsv'
        with summ.open('w', newline='') as fh:
            writer = csv.DictWriter(fh, fieldnames=summary_rows[0].keys(), delimiter='\t')
            writer.writeheader(); [writer.writerow(r) for r in summary_rows]
        log(f"Kraken summary -> {summ}")
        return summ
    return None

# ===================== Post-QC =====================

def locate_clean_fastas() -> Dict[str, Path]:
    """Return mapping of sample -> absolute path to cleaned FASTA."""
    mapping={}
    clean_dir = (KRAKEN_DIR / 'clean').resolve()
    if not clean_dir.exists():
        return mapping
    for fa in clean_dir.glob('*_clean.fasta'):
        sample = fa.name.rsplit('_clean.fasta', 1)[0]
        mapping[sample]=fa.resolve()
    return mapping

def wait_for_clean_fastas(cfg: PipelineConfig) -> Dict[str, Path]:
    """Optionnellement attendre (poll) les fichiers FASTA nettoyés (utile si Kraken exécuté externe / en file d'attente)."""
    deadline = time.time() + cfg.wait_clean_secs
    cleaned = locate_clean_fastas()
    if cleaned or cfg.wait_clean_secs <= 0:
        return cleaned
    log(f'Attente jusqu\'à {cfg.wait_clean_secs}s pour les FASTA nettoyés (sorties Kraken)...')
    while time.time() < deadline:
        time.sleep(5)
        cleaned = locate_clean_fastas()
        if cleaned:
            break
    if not cleaned:
        warn('Délai d\'attente dépassé pour les FASTA nettoyés.')
    return cleaned

def run_post_qc() -> Optional[Path]:
    cleaned=locate_clean_fastas()
    if not cleaned:
        warn('Aucun FASTA nettoyé détecté -> saut de la post-QC')
        return None
    rows=[]
    for sample, path in cleaned.items():
        lengths=[]; total_bp=0
        for cid, seq in fasta_iter(path):
            L=len(seq); lengths.append(L); total_bp+=L
        rows.append({'sample': sample,'n_contigs_clean': len(lengths),'total_bp_clean': total_bp,'N50_clean': n50(lengths)})
    out = REPORTS / 'post_decontam_qc_metrics.tsv'
    with out.open('w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=rows[0].keys(), delimiter='\t')
        writer.writeheader(); [writer.writerow(r) for r in rows]
    log(f"Post-QC metrics -> {out}")
    return out

# ===================== QUAST =====================

def have_quast() -> bool:
    return find_quast() is not None

def _quast_bin():
    return find_quast() or 'quast'

def run_quast(genomes: Dict[str, Path], reference: Optional[Path], threads: int, suffix: str, force: bool=False) -> Optional[Path]:
    if not have_quast():
        warn('QUAST not found')
        return None
    outdir = RESULTS / f'quast_{suffix}'
    outdir.mkdir(parents=True, exist_ok=True)
    transposed = outdir / 'transposed_report.tsv'
    if transposed.exists() and not force:
        log(f'QUAST cache detected ({suffix})')
        return outdir
    
    # Correction: si la référence est un chemin mais pas une clé de génome, l'accepter
    ref_path = None
    if isinstance(reference, Path) and reference.exists():
        ref_path = reference
    elif isinstance(reference, str) and Path(reference).exists():
        ref_path = Path(reference)
    elif reference in genomes:
        ref_path = genomes[reference]

    fasta_files=[str(p) for p in genomes.values() if p.exists()]
    if not fasta_files:
        warn('No FASTA for QUAST')
        return None
    labels=','.join(genomes.keys())
    cmd=[_quast_bin(), '-o', str(outdir), '-t', str(threads), '--labels', labels, '-m', '100']
    if ref_path:
        cmd += ['-r', str(ref_path)]
    cmd += fasta_files
    r=run_cmd(cmd, capture=True)
    if r.returncode!=0:
        err_preview = ''
        if r.stderr:
            try:
                err_preview = ' | '.join(r.stderr.splitlines()[:3])
            except Exception:
                err_preview = r.stderr[:200]
        warn(f"QUAST failed (exit {r.returncode}): {err_preview}")
        return None
    return outdir

def parse_quast(outdir: Path, suffix: str) -> Optional[Path]:
    transposed = outdir / 'transposed_report.tsv'
    if not transposed.exists():
        warn('Missing transposed_report.tsv')
        return None
    import pandas as pd
    try:
        df=pd.read_csv(transposed, sep='\t')
    except Exception as e:
        warn(f'Cannot read QUAST report: {e}')
        return None
    feat_col=df.columns[0]
    idx=df.set_index(feat_col)
    desired={
        '# contigs (>= 0 bp)': f'quast_contigs_{suffix}',
        'Total length (>= 0 bp)': f'quast_total_len_{suffix}',
        'GC (%)': f'quast_GC_{suffix}',
        'N50': f'quast_N50_{suffix}',
        'Largest contig': f'quast_largest_contig_{suffix}',
    }
    rows=[]
    for sample in [c for c in idx.columns]:
        row={'sample': sample}
        for feat, new_name in desired.items():
            if feat in idx.index:
                val=idx.loc[feat, sample]
            else:
                # fallback partial
                partial=[i for i in idx.index if feat.split('(')[0].strip() in i]
                val=idx.loc[partial[0], sample] if partial else None
            row[new_name]=val
        rows.append(row)
    import pandas as _pd
    out = RESULTS / f'quast_summary_{suffix}.tsv'
    _pd.DataFrame(rows).to_csv(out, sep='\t', index=False)
    log(f'QUAST summary -> {out}')
    return out

# ===================== BUSCO =====================

def have_busco():
    return shutil.which('busco') is not None

# Replace the existing run_busco function with this improved version
def run_busco(genomes: Dict[str, Path], lineage: str, threads: int, suffix: str,
              force=False, auto_lineage: bool=False) -> Optional[Path]:
    """
    Exécute BUSCO:
      - Mode dataset explicite (-l lineage) OU auto-lineage (--auto-lineage*)
      - Fallback automatique: si lineage explicite invalide => relance en auto-lineage
    """
    if not have_busco():
        warn('BUSCO not found')
        return None
    outdir = (RESULTS / f'busco_{suffix}').resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    summary = outdir / f'busco_summary_{suffix}.tsv'
    if summary.exists() and not force:
        log(f'BUSCO cache detected ({suffix})')
        return summary

    # Sécuriser variables d'environnement relatives
    for var in ('AUGUSTUS_CONFIG_PATH','BUSCO_CONFIG_FILE','BUSCO_DOWNLOADS'):
        v = os.environ.get(var)
        if v and not v.startswith('/'):
            abs_v = str(Path(v).resolve())
            os.environ[var] = abs_v
            warn(f'{var} converti en absolu: {abs_v}')

    patt = re.compile(r'C:(?P<C>[0-9.]+)%\[S:(?P<S>[0-9.]+)%,D:(?P<D>[0-9.]+)%],F:(?P<F>[0-9.]+)%,M:(?P<M>[0-9.]+)%')
    rows = []

    # Résolution locale du lineage si nom simple + BUSCO_DOWNLOADS
    lineage_path = Path(lineage)
    lineage_is_dir = lineage_path.exists() and lineage_path.is_dir()
    lineage_arg = str(lineage_path.resolve()) if lineage_is_dir else lineage
    if (not lineage_is_dir) and not auto_lineage:
        dl_base = os.environ.get('BUSCO_DOWNLOADS')
        if dl_base:
            cand = Path(dl_base) / 'lineages' / lineage
            if cand.exists():
                lineage_arg = str(cand.resolve())
                lineage_is_dir = True
                log(f'Using local BUSCO lineage: {lineage_arg}')

    def _build_cmd(fasta_abs: Path, tag: str, use_auto: bool, force_flag: bool):
        cmd = [
            'busco',
            '-i', str(fasta_abs),
            '-o', tag,
            '--out_path', str(outdir),
            '-m', 'genome',
            '--cpu', str(threads)
        ]
        if use_auto:
            cmd.append('--auto-lineage')
        else:
            cmd.extend(['-l', lineage_arg])
            if lineage_is_dir:
                cmd.append('--offline')
        if force_flag:
            cmd.append('-f')
        return cmd

    for sample, fasta in genomes.items():
        try:
            fasta_abs = fasta.resolve(strict=True)
        except FileNotFoundError:
            warn(f'BUSCO skip {sample}: FASTA missing -> {fasta}')
            continue
        if not fasta_abs.is_file():
            warn(f'BUSCO skip {sample}: not a file -> {fasta_abs}')
            continue

        sample_tag = f'busco_{sample}_{suffix}'
        sample_out = outdir / sample_tag
        stderr_log = outdir / f'{sample_tag}.stderr.log'

        need_run = True # Toujours vérifier, le cache est géré au début de la fonction
        invalid_lineage_retry = False
        attempted_force_retry = False

        while True:
            # Vérifier si le résumé existe déjà pour ce sample
            if not force and (outdir / f'short_summary.specific.{lineage}.{sample_tag}.txt').exists():
                 log(f'BUSCO reuse existing summary for {sample}')
                 need_run = False

            if need_run:
                mode_auto = auto_lineage or invalid_lineage_retry
                mode_label = 'auto-lineage' if mode_auto else f'lineage={lineage_arg}'
                force_flag = force or attempted_force_retry
                log(f'BUSCO run ({suffix}) sample={sample} ({mode_label}){" [+force]" if force_flag else ""}')
                cmd = _build_cmd(fasta_abs, sample_tag, mode_auto, force_flag)
                r = run_cmd(cmd, capture=True)
                if r.stderr:
                    stderr_log.write_text(r.stderr)
                if r.returncode != 0:
                    stderr_low = (r.stderr or '').lower()
                    # Amélioration du log d'erreur: montre les dernières lignes non vides
                    error_lines = [line for line in (r.stderr or '').splitlines() if line.strip()]
                    preview = ' | '.join(error_lines[-8:]) # Montre les 8 dernières lignes pertinentes
                    if (not mode_auto) and ('not a valid option for \'lineages\'' in stderr_low):
                        warn(f'Lineage invalide pour {sample}, bascule en auto-lineage.')
                        invalid_lineage_retry = True
                        continue
                    if ('already exists' in stderr_low) and not attempted_force_retry:
                        warn(f'Dossier existant pour {sample}, relance avec -f.')
                        attempted_force_retry = True
                        continue
                    warn(f'BUSCO failed {sample} (exit {r.returncode}) preview: {preview}')
                    break
            else:
                # Si on saute le run, on doit quand même parser le résultat existant
                pass

            # Parsing
            short_files = list(sample_out.glob('short_summary*.txt'))
            if not short_files:
                warn(f'No BUSCO short_summary for {sample} in {sample_out}')
                break
            try:
                content = short_files[0].read_text()
            except Exception as e:
                warn(f'Cannot read BUSCO summary for {sample}: {e}')
                break
            m = patt.search(content)
            if not m:
                line_match = [l for l in content.splitlines() if l.strip().startswith('C:')]
                warn(f'BUSCO summary parse fail {sample} preview={line_match[:1]}')
                break
            rows.append({
                'sample': sample,
                f'busco_C_pct_{suffix}': m.group('C'),
                f'busco_S_pct_{suffix}': m.group('S'),
                f'busco_D_pct_{suffix}': m.group('D'),
                f'busco_F_pct_{suffix}': m.group('F'),
                f'busco_M_pct_{suffix}': m.group('M')
            })
            break  # succès -> sortir boucle while
    if not rows:
        warn('BUSCO produced no parsable summaries.')
        return None
    import pandas as pd
    pd.DataFrame(rows).to_csv(summary, sep='\t', index=False)
    log(f'BUSCO summary -> {summary}')
    return summary

# ===================== Fusion QC =====================

def merge_qc(pre: Optional[Path], post: Optional[Path]) -> Optional[Path]:
    """Fusion pré/post avec logique robuste de parsing (reprend analysis_pipeline)."""
    if not pre or not pre.exists():
        return None
    if pd is None:
        warn('pandas indisponible pour fusion QC')
        return None

    def _read_any(p: Path) -> pd.DataFrame:
        for sep in ("\t", ",", None):
            try:
                if sep is None:
                    df = pd.read_csv(p, sep=None, engine='python')
                else:
                    df = pd.read_csv(p, sep=sep)
                if df.empty:
                    continue
                df.columns = [c.strip() for c in df.columns]
                if 'sample' not in df.columns and df.columns:
                    first = df.columns[0]
                    if first.lower().startswith('sample'):
                        df = df.rename(columns={first: 'sample'})
                if 'sample' in df.columns:
                    return df
            except Exception:
                continue
        return pd.DataFrame()

    pre_df = _read_any(pre)
    if pre_df.empty or 'sample' not in pre_df.columns:
        warn("Pré-QC sans colonne sample: abandon fusion")
        return None

    if post and post.exists():
        post_df = _read_any(post)
        if post_df.empty or 'sample' not in post_df.columns:
            warn('Post-QC non lisible: utilisation pré seul')
            merged = pre_df.copy()
        else:
            merged = pre_df.merge(post_df, on='sample', how='left')
            if {'n_contigs','n_contigs_clean'} <= set(merged.columns):
                merged['contigs_reduction_%'] = (1 - (merged['n_contigs_clean']/merged['n_contigs']))*100
            if {'total_bp','total_bp_clean'} <= set(merged.columns):
                merged['total_bp_reduction_%'] = (1 - (merged['total_bp_clean']/merged['total_bp']))*100
            if {'N50','N50_clean'} <= set(merged.columns):
                merged['N50_change_%'] = ((merged['N50_clean'] - merged['N50'])/merged['N50']*100).replace([float('inf'), -float('inf')], 0)
    else:
        merged = pre_df

    out = REPORTS / 'qc_compare_pre_post.tsv'
    merged.to_csv(out, sep='\t', index=False)
    log(f'QC compare -> {out}')
    return out

# ===================== Quality Flags =====================

def compute_quality_flags(pre: Optional[Path], post: Optional[Path], cfg: PipelineConfig) -> Optional[Path]:
    import pandas as pd
    if not pre or not pre.exists(): return None
    pre_df=pd.read_csv(pre, sep='\t')
    post_df=pd.read_csv(post, sep='\t') if post and post.exists() else None
    rows=[]
    for _, row in pre_df.iterrows():
        sample=row['sample']
        if post_df is not None and sample in set(post_df['sample']):
            prow=post_df[post_df['sample']==sample].iloc[0]
            total_bp=prow.get('total_bp_clean', prow.get('total_bp'))
            n_contigs=prow.get('n_contigs_clean', prow.get('n_contigs'))
            N50v=prow.get('N50_clean', prow.get('N50'))
        else:
            total_bp=row.get('total_bp'); n_contigs=row.get('n_contigs'); N50v=row.get('N50')
        gc_val=row.get('GC_mean')
        flags=[]
        size_dev_pct=(abs(total_bp-cfg.expected_size)/cfg.expected_size*100) if total_bp else None
        gc_dev=abs(gc_val-cfg.expected_gc) if gc_val else None
        if size_dev_pct and size_dev_pct>cfg.size_tolerance_pct: flags.append('SIZE_DEVIATION')
        if gc_dev and gc_dev>cfg.gc_tolerance: flags.append('GC_DEVIATION')
        if n_contigs and n_contigs>2: flags.append('FRAGMENTED')
        if not flags: flags=['OK']
        rows.append({
            'sample': sample,
            'final_contigs': n_contigs,
            'final_total_bp': total_bp,
            'final_N50': N50v,
            'final_GC_mean': gc_val,
            'size_deviation_pct': f"{size_dev_pct:.2f}" if size_dev_pct else '',
            'gc_deviation_pct_points': f"{gc_dev:.2f}" if gc_dev else '',
            'quality_flags': ';'.join(flags)
        })
    out=REPORTS/'quality_flags.tsv'
    import pandas as pd
    pd.DataFrame(rows).to_csv(out, sep='\t', index=False)
    log(f'Quality flags -> {out}')
    return out

# ===================== Variants (placeholder) =====================

def normalize_filter_vcfs():
    warn('Variant processing placeholder - integrate bcftools operations here.')
    return None

# ===================== Pangenome Stubs =====================

def create_pansn_combined_fasta(cfg: PipelineConfig, output_path: Path, ref_name: str = 'ref') -> Optional[Path]:
    """Create combined FASTA with PanSN naming: sample#haplotype#contig.
    
    This function:
    1. Adds the reference genome (from external file or samples) with name 'ref#1#genome'
    2. Renames all sample contigs to follow PanSN format: sampleName#1#ctgN
    3. Orders sequences with reference first
    
    Returns the path to the created file, or None on failure.
    """
    log(f'Creating PanSN-formatted combined FASTA -> {output_path}')
    
    # Determine which genomes to use (cleaned or raw)
    cleaned = locate_clean_fastas()
    use_genomes = cleaned if cleaned else cfg.genomes
    if cleaned:
        log(f'Using {len(cleaned)} cleaned genomes for PanSN formatting')
    else:
        log('No cleaned genomes found -> using raw genomes')
    
    # Prepare all genomes including reference
    all_genomes = use_genomes.copy()
    ref_key = None
    
    # Add external reference if provided
    if cfg.pangenome_reference_fasta:
        ref_path = Path(cfg.pangenome_reference_fasta)
        if ref_path.exists():
            ref_key = ref_name
            all_genomes[ref_key] = ref_path
            log(f'Using external reference: {ref_path} (key={ref_key})')
        else:
            warn(f'Reference FASTA not found: {ref_path}')
            return None
    elif cfg.pangenome_reference_key:
        # Use a sample as reference
        if cfg.pangenome_reference_key in all_genomes:
            ref_key = cfg.pangenome_reference_key
            log(f'Using sample as reference: {ref_key}')
        else:
            warn(f'Reference key not found: {cfg.pangenome_reference_key}')
            ref_key = next(iter(all_genomes.keys())) if all_genomes else None
    else:
        # Default to first genome
        ref_key = next(iter(all_genomes.keys())) if all_genomes else None
        if ref_key:
            log(f'No reference specified, using first genome: {ref_key}')
    
    if not ref_key:
        warn('No reference or genomes available')
        return None
    
    # Order: reference first, then other samples
    ordered_keys = [ref_key] + [k for k in all_genomes.keys() if k != ref_key]
    
    # Track contig counts per sample
    contig_counts = {}
    sequences_written = 0
    
    try:
        with output_path.open('w') as fh:
            for sample_key in ordered_keys:
                path = all_genomes[sample_key]
                if not path.exists():
                    warn(f'Genome file not found (skip): {path}')
                    continue
                
                # Special handling for reference: name it 'ref#1#genome'
                if sample_key == ref_key:
                    for cid, seq in fasta_iter(path):
                        new_name = f'{ref_name}#1#genome'
                        fh.write(f'>{new_name}\n{seq}\n')
                        sequences_written += 1
                        log(f'  Reference: {cid} -> {new_name}')
                        break  # Only take first sequence for reference
                else:
                    # Sample contigs: sampleName#1#ctgN
                    if sample_key not in contig_counts:
                        contig_counts[sample_key] = 0
                    
                    for cid, seq in fasta_iter(path):
                        contig_counts[sample_key] += 1
                        new_name = f'{sample_key}#1#ctg{contig_counts[sample_key]}'
                        fh.write(f'>{new_name}\n{seq}\n')
                        sequences_written += 1
                        if contig_counts[sample_key] <= 3:  # Log first few
                            log(f'  {sample_key}: {cid} -> {new_name}')
                    
                    if contig_counts[sample_key] > 3:
                        log(f'  {sample_key}: ... ({contig_counts[sample_key]} total contigs)')
        
        log(f'PanSN FASTA created: {sequences_written} sequences written')
        log(f'Summary: reference + {len(contig_counts)} samples')
        for sample, count in sorted(contig_counts.items()):
            log(f'  {sample}: {count} contig(s)')
        
        return output_path
        
    except Exception as e:
        error(f'Failed to create PanSN FASTA: {e}')
        return None

def run_cactus_pipeline(cfg: PipelineConfig) -> Optional[Path]:
    """Execute Minigraph-Cactus Snakemake pipeline if available.

    Strategy:
      1. Verify presence of cactus-snakemake directory & snakemake binary.
      2. Create minimal config YAML (or reuse provided one) in results/minigraph_cactus.
      3. Run snakemake for main cactus workflow.
      4. Detect final HAL (or other expected artifact) and return its path.
    """
    pipeline_dir = Path(cfg.cactus_pipeline_dir)
    if not pipeline_dir.exists():
        warn(f'Cactus pipeline directory not found: {pipeline_dir}')
        return None
    if not which('snakemake'):
        warn('snakemake not found in PATH - cannot run cactus pipeline')
        return None
    work_root = RESULTS / 'minigraph_cactus'
    work_root.mkdir(parents=True, exist_ok=True)
    config_dir = work_root
    
    # Create PanSN-formatted combined FASTA for cactus
    combined_fa = work_root / 'input_genomes_pansn.fa'
    pansn_result = create_pansn_combined_fasta(cfg, combined_fa, ref_name='ref')
    if not pansn_result:
        warn('Failed to create PanSN-formatted FASTA for cactus')
        return None
    
    # Determine config file
    if cfg.cactus_config_path and Path(cfg.cactus_config_path).exists():
        config_file = Path(cfg.cactus_config_path)
    else:
        config_file = config_dir / 'auto_config.yaml'
        # Build a minimal config using the PanSN-formatted file
        output_dir = work_root / 'output'
        output_dir.mkdir(parents=True, exist_ok=True)
        tmp_dir = work_root / 'tmp'
        tmp_dir.mkdir(parents=True, exist_ok=True)
        seq_lines = [
            '# Auto-generated cactus config with PanSN formatting',
            'genomes:',
            f'  combined: {combined_fa.resolve()}',
            'reference: ref',
            f'output_dir: {output_dir.resolve()}',
            f'tmp_dir: {tmp_dir.resolve()}',
            f'output_prefix: {cfg.cactus_output_prefix}',
            f'threads: {cfg.pangenome_threads}'
        ]
        config_file.write_text('\n'.join(seq_lines) + '\n')
    expected_hal = work_root / 'output' / f"{cfg.cactus_output_prefix}.hal"
    if expected_hal.exists() and cfg.cactus_skip_if_exists:
        log(f'Cactus output HAL already exists -> {expected_hal} (skip)')
        return expected_hal
    snakefile = pipeline_dir / 'cactus.smk'
    if not snakefile.exists():
        warn(f'Snakefile cactus.smk not found in {pipeline_dir}')
        return None
    cmd = [
        'snakemake',
        '-s', str(snakefile),
        '-j', str(cfg.pangenome_threads),
        '--configfile', str(config_file),
        '--directory', str(work_root)
    ]
    start = time.time()
    r = run_cmd(cmd, capture=True)
    if r.returncode != 0:
        warn(f'Cactus pipeline failed (exit {r.returncode})')
        if r.stderr:
            warn(r.stderr.splitlines()[:5])
        return None
    duration = time.time() - start
    log(f'Cactus pipeline completed in {duration/60:.1f} min')
    if expected_hal.exists():
        log(f'Found HAL: {expected_hal}')
        # Optional: convert HAL to basic metrics if halStats exists
        if which('halStats'):
            stats_out = expected_hal.with_suffix('.hal.stats.txt')
            rs = run_cmd(['halStats', str(expected_hal)], capture=True)
            if rs.returncode == 0 and rs.stdout:
                stats_out.write_text(rs.stdout)
                log(f'halStats -> {stats_out}')
        return expected_hal
    else:
        warn('Expected HAL not found after pipeline run')
        return None

def _count_gfa_metrics(gfa_path: Path) -> Dict[str, int]:
    metrics = {'segments':0,'links':0,'paths':0}
    try:
        with gfa_path.open() as fh:
            for line in fh:
                if not line or line.startswith('#'): continue
                t=line[0]
                if t=='S': metrics['segments']+=1
                elif t=='L': metrics['links']+=1
                elif t=='P': metrics['paths']+=1
    except Exception as e:
        warn(f'GFA metric parsing failed: {e}')
    return metrics

def run_pggb_pipeline(cfg: PipelineConfig) -> Optional[Path]:
    """Execute PGGB pangenome graph construction with PanSN naming.

    Steps:
      1. Create PanSN-formatted combined FASTA with reference.
      2. Build PGGB command with user parameters + variant calling.
      3. Run unless GFA already present (caching).
      4. Summarize simple metrics (segments, links, paths).
    """
    if not which('pggb'):
        warn('pggb not found in PATH')
        return None
    pggb_root = RESULTS / 'pggb'
    out_dir = pggb_root / 'output'
    out_dir.mkdir(parents=True, exist_ok=True)
    combined_fa = pggb_root / 'input_genomes_pansn.fa'
    
    # Create PanSN-formatted FASTA (always recreate to ensure correct format)
    pansn_result = create_pansn_combined_fasta(cfg, combined_fa, ref_name='ref')
    if not pansn_result:
        warn('Failed to create PanSN-formatted FASTA')
        return None
    
    # Detect existing GFA
    existing_gfas = list(out_dir.glob('*.gfa'))
    if existing_gfas and cfg.pggb_skip_if_exists:
        log(f'PGGB GFA already present -> {existing_gfas[0]} (skip)')
        return existing_gfas[0]
    
    # Build PGGB command with -V ref:100 for variant calling
    cmd = [
        'pggb', '-i', str(combined_fa), '-o', str(out_dir),
        '-t', str(cfg.pggb_threads), '-s', str(cfg.pggb_segment),
        '-p', str(cfg.pggb_percent_id), '-n', str(cfg.pggb_passes),
        '-V', 'ref:100'  # Enable variant calling with reference 'ref' at 100% coverage
    ]
    if cfg.pggb_extra:
        cmd.extend(cfg.pggb_extra.split())
    start = time.time()
    r = run_cmd(cmd, capture=True)
    if r.returncode != 0:
        warn(f'PGGB failed (exit {r.returncode})')
        if r.stderr:
            # Amélioration du log d'erreur: montre les dernières lignes non vides
            error_lines = [line for line in r.stderr.splitlines() if line.strip()]
            warn("=== PGGB Stderr (last 10 lines) ===")
            for line in error_lines[-10:]:
                warn(line)
            warn("====================================")
        return None
    duration = time.time() - start
    log(f'PGGB completed in {duration/60:.1f} min')
    gfas = list(out_dir.glob('*.gfa'))
    if not gfas:
        warn('No GFA produced by PGGB')
        return None
    gfa_main = gfas[0]
    metrics = _count_gfa_metrics(gfa_main)
    metrics_path = out_dir / 'graph_metrics.tsv'
    with metrics_path.open('w') as fh:
        fh.write('metric\tvalue\n')
        for k,v in metrics.items():
            fh.write(f'{k}\t{v}\n')
    log(f'GFA metrics -> {metrics_path}')
    return gfa_main

# ===================== Advanced Analyses (stubs) =====================

def run_phylogeny(): warn('Phylogeny stub - requires alignment + tree builder.')

def run_synteny(): warn('Synteny stub - requires alignment tool like minimap/mummer.')

def run_functional_analysis(): warn('Functional impact stub - requires variant annotation.')

# ===================== Global Summary =====================

def build_global_summary():
    import pandas as pd
    dfs={}
    for name in ['pre_decontam_qc_metrics.tsv','post_decontam_qc_metrics.tsv','qc_compare_pre_post.tsv','quality_flags.tsv']:
        p=REPORTS/name
        if p.exists(): dfs[name]=pd.read_csv(p, sep='\t')
    if not dfs:
        warn('No tables to summarize.')
        return None
    base=None
    for key in ['qc_compare_pre_post.tsv','pre_decontam_qc_metrics.tsv']:
        if key in dfs: base=dfs[key]; break
    merged=base
    for k, d in dfs.items():
        if d is merged: continue
        if 'sample' in d.columns and 'sample' in merged.columns:
            merged=merged.merge(d, on='sample', how='left', suffixes=('','_'+k.replace('.tsv','')))
    out=REPORTS/'global_summary.tsv'
    merged.to_csv(out, sep='\t', index=False)
    log(f'Global summary -> {out}')
    return out

# ===================== Orchestrateur =====================
STEP_FUNCS: Dict[str, Callable[[PipelineConfig], Optional[Path]]] = {}

def register_step(name: str):
    def decorator(func):
        STEP_FUNCS[name]=func
        return func
    return decorator

@register_step('pre-qc')
def step_pre_qc(cfg: PipelineConfig): return run_pre_qc(cfg)

@register_step('kraken')
def step_kraken(cfg: PipelineConfig): return submit_kraken_jobs(cfg)

@register_step('post-qc')
def step_post_qc(cfg: PipelineConfig): return run_post_qc()

@register_step('quast-pre')
def step_quast_pre(cfg: PipelineConfig):
    ref_key_or_path = cfg.reference_key or cfg.pangenome_reference_key
    ref = None
    if ref_key_or_path:
        p = Path(ref_key_or_path)
        if p.is_file():
            ref = p
        elif ref_key_or_path in cfg.genomes:
            ref = cfg.genomes[ref_key_or_path]
    if not ref:
        ref = next(iter(cfg.genomes.values()))

    outdir = run_quast(cfg.genomes, ref, cfg.quast_threads, 'pre', cfg.force_quast)
    return parse_quast(outdir, 'pre') if outdir else None

@register_step('busco-pre')
def step_busco_pre(cfg: PipelineConfig):
    return run_busco(cfg.genomes, cfg.busco_lineage, cfg.busco_threads, 'pre', cfg.force_busco, cfg.busco_auto_lineage)

@register_step('quast-post')
def step_quast_post(cfg: PipelineConfig):
    cleaned = wait_for_clean_fastas(cfg)
    if not cleaned:
        warn('No cleaned FASTA for QUAST post')
        return None
    
    ref_key_or_path = cfg.reference_key or cfg.pangenome_reference_key
    ref = None
    if ref_key_or_path:
        p = Path(ref_key_or_path)
        if p.is_file():
            ref = p
        elif ref_key_or_path in cfg.genomes:
            ref = cfg.genomes[ref_key_or_path]
    if not ref:
        ref = next(iter(cfg.genomes.values()))

    # Ensure absolute
    cleaned_abs = {k: v.resolve() for k,v in cleaned.items()}
    outdir = run_quast(cleaned_abs, ref.resolve() if ref else None, cfg.quast_threads, 'post', cfg.force_quast)
    return parse_quast(outdir, 'post') if outdir else None

@register_step('busco-post')
def step_busco_post(cfg: PipelineConfig):
    cleaned = wait_for_clean_fastas(cfg)
    if not cleaned:
        warn('No cleaned FASTA for BUSCO post')
        return None
    cleaned_abs = {k: v.resolve() for k,v in cleaned.items()}
    return run_busco(cleaned_abs, cfg.busco_lineage, cfg.busco_threads, 'post', cfg.force_busco, cfg.busco_auto_lineage)

@register_step('qc-merge')
def step_qc_merge(cfg: PipelineConfig):
    pre = REPORTS/'pre_decontam_qc_metrics.tsv'
    post = REPORTS/'post_decontam_qc_metrics.tsv'
    return merge_qc(pre, post)

@register_step('quality-flags')
def step_quality_flags(cfg: PipelineConfig):
    pre=REPORTS/'pre_decontam_qc_metrics.tsv'; post=REPORTS/'post_decontam_qc_metrics.tsv'
    return compute_quality_flags(pre, post, cfg)

@register_step('variants')
def step_variants(cfg: PipelineConfig): return normalize_filter_vcfs()

@register_step('pangenome-cactus')
def step_cactus(cfg: PipelineConfig): return run_cactus_pipeline(cfg)

@register_step('pangenome-pggb')
def step_pggb(cfg: PipelineConfig): return run_pggb_pipeline(cfg)

@register_step('phylogeny')
def step_phylogeny(cfg: PipelineConfig): run_phylogeny()

@register_step('synteny')
def step_synteny(cfg: PipelineConfig): run_synteny()

@register_step('functional')
def step_functional(cfg: PipelineConfig): run_functional_analysis()

@register_step('summary')
def step_summary(cfg: PipelineConfig): return build_global_summary()

ALL_ORDER = [
    'pre-qc','kraken','post-qc','quast-pre','busco-pre','quast-post','busco-post',
    'qc-merge','quality-flags','variants','pangenome-cactus','pangenome-pggb',
    'phylogeny','synteny','functional','summary'
]

# ===================== CLI =====================

def parse_args():
    ap=argparse.ArgumentParser(description='Unified KHV complete pipeline')
    ap.add_argument('--genomes', nargs='+', help='sample:path ...')
    ap.add_argument('--reference', help='Reference sample key')
    ap.add_argument('--steps', help='Comma-separated steps (see --list-steps)')
    ap.add_argument('--all', action='store_true', help='Run all steps')
    ap.add_argument('--kraken-db', default=DB_DEFAULT)
    ap.add_argument('--target-taxids', help='Comma-separated taxid(s) to retain after Kraken (optional)')
    ap.add_argument('--keep-unclassified', action='store_true', help='Retain unclassified contigs in cleaned FASTA')
    ap.add_argument('--busco-lineage', default='viruses_odb10')
    ap.add_argument('--busco-auto-lineage', action='store_true',
                    help='Utiliser BUSCO --auto-lineage (ignore --busco-lineage si présent)')
    ap.add_argument('--quast-threads', type=int, default=int(os.environ.get('SLURM_CPUS_PER_TASK','4')))
    ap.add_argument('--busco-threads', type=int, default=int(os.environ.get('SLURM_CPUS_PER_TASK','4')))
    ap.add_argument('--force-quast', action='store_true')
    ap.add_argument('--force-busco', action='store_true')
    ap.add_argument('--no-plots', action='store_true')
    ap.add_argument('--skip-kraken-if-present', action='store_true')
    ap.add_argument('--expected-size', type=int, default=EXPECTED_DEFAULT_SIZE)
    ap.add_argument('--expected-gc', type=float, default=EXPECTED_DEFAULT_GC)
    ap.add_argument('--size-tolerance-pct', type=float, default=12.0)
    ap.add_argument('--gc-tolerance', type=float, default=5.0)
    # Pangenome specific
    ap.add_argument('--pangenome-threads', type=int, default=int(os.environ.get('SLURM_CPUS_PER_TASK','8')))
    ap.add_argument('--cactus-pipeline-dir', default='cactus-snakemake')
    ap.add_argument('--cactus-output-prefix', default='khv')
    ap.add_argument('--cactus-config', help='Path to existing cactus config file')
    ap.add_argument('--no-cactus-skip', action='store_true', help='Do not skip cactus if HAL exists')
    ap.add_argument('--pggb-segment', type=int, default=5000)
    ap.add_argument('--pggb-pid', type=int, default=70)
    ap.add_argument('--pggb-passes', type=int, default=10)
    ap.add_argument('--pggb-threads', type=int, default=int(os.environ.get('SLURM_CPUS_PER_TASK','8')))
    ap.add_argument('--pggb-extra', default='', help='Extra raw flags appended to pggb command')
    ap.add_argument('--no-pggb-skip', action='store_true', help='Do not skip if GFA exists')
    ap.add_argument('--pangenome-reference', help='Sample key to use systematically as pangenome reference (ordering + annotation)')
    ap.add_argument('--pangenome-reference-fasta', help='Path to external reference FASTA file for pangenome (e.g., clean_fasta/ref/KHV-U_trunc.fasta)')
    ap.add_argument('--list-steps', action='store_true')
    ap.add_argument('--wait-clean-secs', type=int, default=0,
                    help='Max seconds to wait/poll for Kraken cleaned FASTA before post analyses')
    return ap.parse_args()

def main():
    args=parse_args()
    if args.list_steps:
        print('Available steps (execution order):')
        for s in ALL_ORDER: print(' -', s)
        sys.exit(0)
    
    if not args.genomes:
        error('--genomes is required'); sys.exit(1)
    
    enrich_path()
    genomes={}
    for entry in args.genomes:
        if ':' not in entry:
            error(f'Invalid genome spec: {entry}'); sys.exit(1)
        k,p = entry.split(':',1)
        genomes[k]=Path(p)
    cfg=PipelineConfig(
        genomes=genomes,
        reference_key=args.reference,
        pangenome_reference_key=args.pangenome_reference or args.reference,
        pangenome_reference_fasta=args.pangenome_reference_fasta,
        kraken_db=args.kraken_db,
        target_taxids=[int(x) for x in args.target_taxids.split(',')] if args.target_taxids else [],
        keep_unclassified=args.keep_unclassified,
        busco_lineage=args.busco_lineage,
        busco_auto_lineage=args.busco_auto_lineage,  # <--- AJOUT
        quast_threads=args.quast_threads,
        busco_threads=args.busco_threads,
        force_quast=args.force_quast,
        force_busco=args.force_busco,
        skip_kraken_if_present=args.skip_kraken_if_present,
        plots=not args.no_plots,
        expected_size=args.expected_size,
        expected_gc=args.expected_gc,
        size_tolerance_pct=args.size_tolerance_pct,
        gc_tolerance=args.gc_tolerance,
        wait_clean_secs=args.wait_clean_secs,
        pangenome_threads=args.pangenome_threads,
        cactus_pipeline_dir=args.cactus_pipeline_dir,
        cactus_output_prefix=args.cactus_output_prefix,
        cactus_config_path=args.cactus_config,
        cactus_skip_if_exists=not args.no_cactus_skip,
        pggb_segment=args.pggb_segment,
        pggb_percent_id=args.pggb_pid,
        pggb_passes=args.pggb_passes,
        pggb_threads=args.pggb_threads,
        pggb_extra=args.pggb_extra,
        pggb_skip_if_exists=not args.no_pggb_skip,
    )
    ensure_dirs()
    check_dependencies()

    if args.all:
        steps=ALL_ORDER.copy()
    elif args.steps:
        steps=[s.strip() for s in args.steps.split(',') if s.strip()]
    else:
        error('Specify --all or --steps')
        sys.exit(1)

    # Dependency enforcement for post steps (quast-post/busco-post need cleaned FASTA)
    need_post = any(s in steps for s in ('quast-post','busco-post'))
    if need_post:
        cleaned = locate_clean_fastas()
        if not cleaned:
            if 'kraken' not in steps:
                # Insert kraken before the first post step (prefer after pre-qc if present)
                insert_idx = 0
                if 'pre-qc' in steps:
                    insert_idx = steps.index('pre-qc') + 1
                steps.insert(insert_idx, 'kraken')
                warn("Added missing dependency step 'kraken' before post analyses.")
            if 'post-qc' not in steps:
                # Place post-qc right after kraken
                k_index = steps.index('kraken')
                steps.insert(k_index + 1, 'post-qc')
                warn("Added missing dependency step 'post-qc' before post analyses.")

    # Deduplicate while preserving order
    seen=set(); ordered=[]
    for s in steps:
        if s not in seen:
            ordered.append(s); seen.add(s)
    steps=ordered

    for step in steps:
        if step not in STEP_FUNCS:
            warn(f'Unknown step: {step} (skip)')
            continue
        log(f"=== Executing step: {step} ===")
        try:
            STEP_FUNCS[step](cfg)
        except KeyboardInterrupt:
            error('Interrupted by user')
            break
        except Exception as e:
            error(f'Step {step} failed: {e}')
    log('Pipeline complete.')

if __name__=='__main__':
    main()
