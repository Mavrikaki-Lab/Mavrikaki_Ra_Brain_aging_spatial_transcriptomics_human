"""Train reference signatures and map spatial cell-type abundances from prepared counts.

CPU workflow; optional comparison against an existing abundance result.
See README.md for input preparation requirements and provenance.
"""
import argparse
import datetime
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import re
import sys
import time


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--spatial', type=Path, required=True, help='Prepared spatial count AnnData')
    parser.add_argument('--reference', type=Path, required=True, help='Prepared single-cell reference count AnnData')
    parser.add_argument('--output-dir', type=Path, required=True, help='New directory; must not exist')
    parser.add_argument('--baseline', type=Path, help='Optional existing spatial abundance result')
    parser.add_argument('--study-input-checks', action='store_true', help='Require the exact study dimensions and category counts')
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--reference-epochs', type=int, default=400)
    parser.add_argument('--mapping-epochs', type=int, default=3000)
    parser.add_argument('--posterior-samples', type=int, default=1000)
    parser.add_argument('--n-cells-per-location', type=float, default=5)
    parser.add_argument('--detection-alpha', type=float, default=20)
    parser.add_argument('--reference-labels', default='celltypes_final')
    parser.add_argument('--spatial-batch', default='sample')
    parser.add_argument('--check-only', action='store_true')
    args = parser.parse_args()
    if min(args.reference_epochs, args.mapping_epochs, args.posterior_samples) < 1:
        parser.error('Epochs and posterior samples must be positive')
    if not (0 < args.n_cells_per_location < float('inf') and 0 < args.detection_alpha < float('inf')):
        parser.error('Cell count and detection priors must be finite and positive')
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=False)
    # Select CPU before importing torch/scvi. This preserves the validated execution mode.
    os.environ['CUDA_VISIBLE_DEVICES'] = ''
    os.environ.setdefault('MPLBACKEND', 'Agg')
    versions = {'python': sys.version, 'executable': sys.executable}
    for package in ['cell2location','scvi-tools','torch','pyro-ppl','anndata','scanpy','numpy','pandas','scipy','lightning','pytorch-lightning']:
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = 'not installed'
    (out/'versions.json').write_text(json.dumps(versions, indent=2))
    print('Output:', out, flush=True)
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp
    import scanpy as sc
    import torch
    import matplotlib.pyplot as plt
    from cell2location.models import RegressionModel, Cell2location
    from scvi import settings
    if torch.cuda.is_available():
        raise RuntimeError('CPU execution was requested but CUDA is still visible')
    settings.seed = args.seed
    config = dict(seed=args.seed,reference_epochs=args.reference_epochs,mapping_epochs=args.mapping_epochs,posterior_samples=args.posterior_samples,
                  reference_labels=args.reference_labels,reference_batch=None,spatial_batch=args.spatial_batch,
                  N_cells_per_location=args.n_cells_per_location,detection_alpha=args.detection_alpha,device='cpu',
                  reference_posterior_batch='all_reference_cells',mapping_posterior_batch=4096)
    (out/'config.json').write_text(json.dumps(config,indent=2))
    paths = {'spatial': args.spatial.resolve(), 'reference': args.reference.resolve()}
    if args.baseline:
        paths['baseline'] = args.baseline.resolve()
    (out/'run_provenance.json').write_text(json.dumps({
        'started_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'command': sys.argv, 'script_sha256': sha256(Path(__file__)),
        'study_input_checks': args.study_input_checks,
        'threads': {k: os.environ.get(k) for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS']}
    }, indent=2))
    for path in paths.values():
        if not path.is_file():
            raise FileNotFoundError(path)
    (out/'inputs.json').write_text(json.dumps({k:{'path':str(p),'sha256':sha256(p)} for k,p in paths.items()},indent=2))
    print('Loading prepared inputs...',flush=True)
    vis = sc.read_h5ad(paths['spatial'])
    ref = sc.read_h5ad(paths['reference'])
    if args.study_input_checks and (vis.shape != (167613,18085) or ref.shape != (287818,1622)):
        raise ValueError(f'Input shapes differ from study inputs: vis={vis.shape}, ref={ref.shape}')
    for ad,col in [(vis,config['spatial_batch']),(ref,config['reference_labels'])]:
        if col not in ad.obs or ad.obs[col].isna().any():
            raise ValueError('Missing annotation: '+col)
    if args.study_input_checks and (vis.obs[config['spatial_batch']].nunique()!=37 or ref.obs[config['reference_labels']].nunique()!=21):
        raise ValueError('Expected 37 spatial samples and 21 reference cell types')
    # Check nonnegative integer-like X in blocks; do not silently select another layer.
    # Prepared counts must be in X, not normalized expression or a fallback layer.
    for tag,ad in [('spatial',vis),('reference',ref)]:
        for start in range(0,ad.n_obs,2048):
            block=ad.X[start:start+2048]
            vals=block.data if sp.issparse(block) else np.asarray(block).ravel()
            if not (np.isfinite(vals).all() and (vals>=0).all() and (np.abs(vals-np.round(vals))<=1e-6).all()):
                raise ValueError(tag+' X is not nonnegative finite integer-like counts')
    ref.var_names=[re.sub(r'\.\d+$','',str(g)) for g in ref.var_names]
    if 'gene_id' not in vis.var:
        raise ValueError('Missing spatial var[gene_id]')
    if vis.var['gene_id'].isna().any():
        raise ValueError('Missing spatial gene IDs')
    vis.var['ENSEMBL']=[re.sub(r'\.\d+$','',str(g)) for g in vis.var['gene_id']]
    vis.var_names=vis.var['ENSEMBL'].astype(object)
    if vis.var_names.hasnans or not vis.var_names.is_unique:
        raise ValueError('Spatial gene IDs must be unique after version stripping; resolve duplicates during preparation')
    if not ref.var_names.is_unique or ref.var_names.hasnans:
        raise ValueError('Duplicate reference genes after version stripping')
    inter=np.intersect1d(ref.var_names,vis.var_names)
    if len(inter)==0:
        raise ValueError("Reference and spatial inputs have no common genes")
    if args.study_input_checks and len(inter)!=1612:
        raise ValueError(f'Expected 1612 common genes, found {len(inter)}')
    np.savetxt(out/'aligned_genes.txt',inter,fmt='%s')
    vis_ids = pd.Index(vis.obs['key'] if 'key' in vis.obs else vis.obs_names)
    if vis_ids.hasnans or vis_ids.astype(str).has_duplicates or (vis_ids.astype(str)=='').any():
        raise ValueError('Missing or ambiguous spatial spot identifiers')
    if min(vis.shape + ref.shape)==0:
        raise ValueError('Empty input matrix')
    labels = set(ref.obs[config['reference_labels']].astype(str))
    if 'baseline' in paths:
        from compare_abundances import extract, compare
        if config['spatial_batch'] != 'sample':
            raise ValueError('Optional comparison currently requires --spatial-batch sample')
        old=sc.read_h5ad(paths['baseline'],backed='r')
        try:
            q=extract(old)
            if set(q.index)!=set(vis_ids.astype(str)) or set(q.columns)!=labels:
                raise ValueError('Baseline spot IDs or cell-type labels differ from the inputs')
            if 'sample' not in old.obs or old.obs['sample'].isna().any():
                raise ValueError('Missing baseline sample annotation')
            samples=pd.Series(old.obs['sample'].astype(str).to_numpy(), index=q.index)
            current=pd.Series(vis.obs['sample'].astype(str).to_numpy(), index=vis_ids.astype(str))
            if not samples.equals(current.loc[samples.index]):
                raise ValueError('Baseline and input sample identities differ')
        finally:
            old.file.close()
    print(f'Input checks passed: {vis.n_obs} spots, {ref.n_obs} reference cells, '
          f'{len(labels)} cell types, {len(inter)} common genes.',flush=True)
    if args.check_only:
        print('PREFLIGHT COMPLETE. No model trained.',flush=True)
        return
    ref=ref[:,inter].copy();vis=vis[:,inter].copy()
    times={}
    def timed(name,fn):
        t=time.monotonic();print('START:',name,flush=True)
        result=fn();times[name]=time.monotonic()-t
        (out/'timings_seconds.json').write_text(json.dumps(times,indent=2))
        print('DONE:',name,round(times[name],1),'seconds',flush=True)
        return result
    RegressionModel.setup_anndata(ref,labels_key=config['reference_labels'],batch_key=config['reference_batch'])
    reg=RegressionModel(ref)
    timed('reference_training',lambda:reg.train(max_epochs=config['reference_epochs']))
    reg.save(str(out/'reference_model'),overwrite=False)
    timed('reference_posterior',lambda:reg.export_posterior(ref,sample_kwargs={'num_samples':config['posterior_samples'],'batch_size':ref.n_obs}))
    if 'means_per_cluster_mu_fg' in ref.varm:
        signatures=ref.varm['means_per_cluster_mu_fg']
    else:
        signatures=ref.uns['mod']['means_per_cluster_mu_fg']
    if not isinstance(signatures,pd.DataFrame):
        signatures=pd.DataFrame(signatures,index=ref.var_names,columns=ref.uns['mod']['factor_names'])
    signatures=signatures.loc[:,signatures.var(0)>0]
    if signatures.shape!=(len(inter),len(labels)):
        raise ValueError(f'Unexpected signature shape: {signatures.shape}')
    signatures.to_csv(out/'reference_signatures.csv')
    # Retain the signature naming and gene order used by the study workflow.
    Cell2location.setup_anndata(vis,batch_key=config['spatial_batch'])
    model=Cell2location(vis,cell_state_df=signatures,N_cells_per_location=config['N_cells_per_location'],detection_alpha=config['detection_alpha'])
    timed('mapping_training',lambda:model.train(max_epochs=config['mapping_epochs'],batch_size=None,train_size=1.0))
    model.save(str(out/'mapping_model'),overwrite=False)
    try:
        model.plot_history(30);plt.savefig(out/'mapping_training_history.png',dpi=150,bbox_inches='tight');plt.close()
    except Exception as e:
        print('Optional history plot failed:',e,flush=True)
    timed('mapping_posterior',lambda:model.export_posterior(vis,sample_kwargs={'num_samples':config['posterior_samples'],'batch_size':min(vis.n_obs,config['mapping_posterior_batch'])}))
    # Posterior failures propagate; never report success without abundance estimates.
    if 'q05_cell_abundance_w_sf' not in vis.obsm:
        raise ValueError('Posterior export did not produce q05 abundances')
    q05 = vis.obsm['q05_cell_abundance_w_sf']
    if np.asarray(q05).shape != (vis.n_obs,len(labels)) or not np.isfinite(np.asarray(q05)).all() or (np.asarray(q05)<0).any():
        raise ValueError('Invalid posterior abundance estimates')
    timed('save_result',lambda:vis.write_h5ad(out/'adata_vis_cell2location.h5ad'))
    if 'baseline' in paths:
        timed('compare_baseline',lambda:compare(vis,paths['baseline'],out/'comparison'))
    try:
        model.plot_QC();plt.savefig(out/'mapping_QC.png',dpi=150,bbox_inches='tight');plt.close()
    except Exception as e:
        print('Optional QC plot failed:',e,flush=True)
    print('CELL2LOCATION COMPLETE\nOutput:',out,flush=True)


if __name__=='__main__':
    main()
