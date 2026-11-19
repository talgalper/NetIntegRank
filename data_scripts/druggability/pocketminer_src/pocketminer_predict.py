#!/usr/bin/env python

### PocketMiner predictions for a batch of structures (called by pocketminer_master.sh) ###
#
# PocketMiner is a trained neural network. It reads the backbone of a protein structure and gives
# every residue a score between 0 and 1: how likely that residue is to take part in a "cryptic"
# pocket, one that is closed in the structure but can open as the protein moves.
#
# usage: pocketminer_predict.py versions
#        pocketminer_predict.py todo <alphafold|swissmodel> <structures_dir> <results_dir>
#        pocketminer_predict.py run  <alphafold|swissmodel> <results_dir> <model_prefix> <pdb> [<pdb> ...]

import os
import sys

## where the authors' PocketMiner code sits (copied unchanged, see upstream_provenance.txt)
script_dir = os.path.dirname(os.path.abspath(__file__))
upstream_src = os.path.join(script_dir, 'upstream', 'src')


## output name for one structure:
## AlphaFold files are named by protein and fragment, SWISS-MODEL files by protein
def output_name(pdb_path, source):
    name = os.path.basename(pdb_path)[:-len('.pdb')]
    if source == 'alphafold':
        parts = name.split('-')                      # AF-P12345-F1-model_v6 -> P12345-F1
        if len(parts) < 3:
            sys.exit(f'{name}.pdb is not named like an AlphaFold file (AF-<protein>-<fragment>-model_<version>.pdb)')
        return parts[1] + '-' + parts[2]
    return name                                      # P12345 -> P12345


## where the per-residue predictions for one structure are saved
def prediction_file(pdb_path, source, results_dir):
    return os.path.join(results_dir, 'txt_results', output_name(pdb_path, source) + '-predictions.txt')


## load TensorFlow and the PocketMiner code
## each worker uses one CPU core, because the master script runs many workers side by side
def load_pocketminer():
    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'         # hide TensorFlow start-up messages
    os.environ['CUDA_VISIBLE_DEVICES'] = ''          # CPU only, so parallel workers never compete for a GPU
    os.environ['OMP_NUM_THREADS'] = '1'
    import warnings
    warnings.filterwarnings('ignore', category=FutureWarning)
    warnings.filterwarnings('ignore', message='Unlikely unit cell vectors')   # mdtraj note on the dummy CRYST1 line, harmless
    import numpy as np
    import tensorflow as tf
    import mdtraj as md
    tf.get_logger().setLevel('ERROR')
    tf.config.threading.set_intra_op_parallelism_threads(1)
    tf.config.threading.set_inter_op_parallelism_threads(1)
    sys.path.insert(0, upstream_src)
    from models import MQAModel
    from validate_performance_on_xtals import process_strucs
    return np, tf, md, MQAModel, process_strucs


## check that the trained weights really are in the model
## method: read the saved weights straight from the checkpoint file and compare them with what the
## model now holds. A wrong checkpoint path does not raise an error in TensorFlow: the model just
## keeps random starting weights and its output is noise, so the run stops here if a single weight differs.
def check_weights_loaded(np, tf, model, model_prefix):
    try:
        reader = tf.train.load_checkpoint(model_prefix)
        saved = [reader.get_tensor(key) for key in reader.get_variable_to_shape_map()
                 if key.startswith('model/') and key.endswith('.ATTRIBUTES/VARIABLE_VALUE')
                 and '.OPTIMIZER_SLOT' not in key]
    except Exception as error:                       # a checkpoint that cannot be read counts as not loaded
        sys.stderr.write(f'FATAL: PocketMiner checkpoint {model_prefix} could not be read: {error}\n')
        sys.exit(255)
    n_weights = len(model.variables)
    n_loaded = sum(any(s.shape == v.shape and np.array_equal(s, v.numpy()) for s in saved)
                   for v in model.variables)
    if n_weights == 0 or n_loaded != n_weights or len(saved) != n_weights:
        sys.stderr.write(f'FATAL: PocketMiner weights not loaded from {model_prefix}: {n_loaded} of {n_weights} '
                         f'model weights match the {len(saved)} saved in the checkpoint\n')
        sys.exit(255)                                # 255 tells xargs to stop the whole run


## mode "versions": print the software versions for run_info.txt
def print_versions():
    np, tf, md, _, _ = load_pocketminer()
    print(f'python: {sys.version.split()[0]}')
    print(f'tensorflow: {tf.__version__}')
    print(f'numpy: {np.__version__}')
    print(f'mdtraj: {md.version.version}')


## mode "todo": list the structures that still need a prediction, one path per line
## structures that already have an output file are left out, which is what lets a stopped run resume
def print_todo(source, structures_dir, results_dir):
    pdbs = sorted(f for f in os.listdir(structures_dir) if f.endswith('.pdb'))
    names = [output_name(f, source) for f in pdbs]
    if len(set(names)) != len(names):
        sys.exit(f'{structures_dir} holds more than one file for the same protein fragment, the outputs would overwrite each other')
    for f in pdbs:
        path = os.path.join(structures_dir, f)
        if not os.path.exists(prediction_file(path, source, results_dir)):
            print(path)


## mode "run": predict every structure in the batch
def run_batch(source, results_dir, model_prefix, pdb_paths):
    np, tf, md, MQAModel, process_strucs = load_pocketminer()

    ## build the network with the settings the authors trained it with (from their xtal_predict.py),
    ## then load the trained weights. The model is built once and reused for the whole batch.
    model = MQAModel(node_features=(8, 50), edge_features=(1, 32), hidden_dim=(16, 100),
                     num_layers=4, dropout=0.1)
    checkpoint = tf.train.Checkpoint(optimizer=tf.keras.optimizers.Adam(), model=model)
    status = checkpoint.restore(model_prefix)
    weights_checked = False

    for pdb_path in pdb_paths:
        out_file = prediction_file(pdb_path, source, results_dir)
        if os.path.exists(out_file):                 # finished in an earlier run
            continue
        try:
            ## read the structure and keep the backbone atoms (N, CA, C, O) of each residue,
            ## then run the network: one score per residue
            X, S, mask = process_strucs([md.load(pdb_path)])
            predictions = np.asarray(model(X, S, mask, train=False, res_level=True))

            ## the weights only exist once the network has run, so they are checked after the first structure
            if not weights_checked:
                check_weights_loaded(np, tf, model, model_prefix)
                weights_checked = True

            ## save one score per line (4 significant figures). The file is written under
            ## a temporary name and renamed when complete, so a half-written file never counts as done.
            np.savetxt(out_file + '.tmp', predictions, fmt='%.4g', delimiter='\n')
            os.replace(out_file + '.tmp', out_file)
        except Exception as error:
            ## a structure that cannot be read or scored is recorded and the batch carries on
            reason = f'{type(error).__name__}: {error}'.replace('\n', ' ')[:300]
            with open(os.path.join(results_dir, 'pocketminer_failed.txt'), 'a') as failed:
                failed.write(f'{output_name(pdb_path, source)}\t{reason}\n')

    del status                                       # keep TensorFlow's end-of-run notes out of the log


## pick the mode from the first argument
if __name__ == '__main__':
    mode = sys.argv[1] if len(sys.argv) > 1 else ''
    if mode == 'versions':
        print_versions()
    elif mode == 'todo' and len(sys.argv) == 5:
        print_todo(sys.argv[2], sys.argv[3], sys.argv[4])
    elif mode == 'run' and len(sys.argv) >= 6:
        run_batch(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5:])
    else:
        sys.exit(__doc__ or 'usage: pocketminer_predict.py versions | todo <source> <structures_dir> <results_dir> '
                            '| run <source> <results_dir> <model_prefix> <pdb> [<pdb> ...]')
