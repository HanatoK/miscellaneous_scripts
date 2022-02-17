#!/usr/bin/env python3
import pandas as pd
import numpy as np
import tensorflow as tf
from tensorflow.keras.models import model_from_json
from sklearn.metrics import explained_variance_score
import matplotlib
matplotlib.use("pgf")
import matplotlib.pyplot as plt
from matplotlib.figure import figaspect
from matplotlib.ticker import AutoMinorLocator


def read_colvars_traj(colvars_traj_filename):
    '''read colvars trajectory of a NANMA simulation with weights'''
    data = pd.read_csv(colvars_traj_filename, delimiter='\s+', comment='#', header=None)
    dihedral_columns = [3, 4]
    dihedral_data = data[dihedral_columns].apply(np.radians)
    sin_phi_psi = dihedral_data.apply(np.sin)
    cos_phi_psi = dihedral_data.apply(np.cos)
    # return a concatenated dataframe
    results = pd.concat([sin_phi_psi, cos_phi_psi, data[5]], axis=1)
    # renumber the column indexes
    for col in range(0, results.shape[1]):
        results.columns.values[col] = str(col)
    return results


def load_model(model_filename, weight_filename):
    '''load tensorflow model from saved files'''
    with open(model_filename, 'r') as json_file:
        model = model_from_json(json_file.read())
        model.load_weights(weight_filename)
        return model


def FVE(model, data):
    input_data = data[[0,1,2,3]]
    output_data = model.predict(input_data, batch_size=input_data.shape[0])
    weights = None
    if data.shape[1] == 5:
        weights = data[4]
    fve = explained_variance_score(y_true=input_data, y_pred=output_data, sample_weight=weights)
    return fve


def plot_fve(models_fve, output_filename):
    plt.rcParams.update({
        "pgf.texsystem": "lualatex",
        "font.family": "serif",  # use serif/main font for text elements
        "text.usetex": True,     # use inline math for ticks
        "pgf.rcfonts": False,    # don't setup fonts from rc parameters
        "axes.labelsize": 26,
        "axes.linewidth": 2.0,
        "font.size": 22,
        "axes.unicode_minus": False,
        "pgf.preamble": '\n'.join([
            "\\usepackage{units}",
            "\\usepackage{metalogo}",
            "\\usepackage{unicode-math}",
            r"\setmathfont{MathJax_Math}",
            r"\setmainfont{FreeSans}",
        ])
    })
    w, h = figaspect(1/1)
    plt.figure(figsize=(w, h))
    colors = matplotlib.cm.get_cmap('Set1')
    num_models = len(models_fve[0])
    x = np.arange(1, num_models + 1)
    for i, y in enumerate(models_fve):
        plt.scatter(x, y, color=colors(i), label=f'Traj {i+1}', alpha=0.7)
    ax = plt.gca()
    ax.set_xlim(0, num_models + 1)
    ax.set_ylim(0, 1)
    ax.tick_params(direction='in', which='major', length=6.0,
                   width=1.0, top=True, right=True, pad=8.0)
    ax.tick_params(direction='in', which='minor', length=3.0,
                   width=1.0, top=True, right=True, pad=8.0)
    ax.xaxis.get_major_formatter()._usetex = False
    ax.yaxis.get_major_formatter()._usetex = False
    ax.xaxis.set_major_locator(plt.MultipleLocator(1))
    ax.yaxis.set_major_locator(plt.MultipleLocator(0.2))
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    plt.xlabel('Round')
    plt.ylabel('Explained variation')
    plt.legend(fontsize=18, handlelength=0.5, handletextpad=0.3,
               handleheight=0.35, labelspacing=0.5, fancybox=False,
               frameon=False, loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig(output_filename, dpi=400, bbox_inches='tight', transparent=False)


if __name__ == '__main__':
    # tf.compat.v1.disable_eager_execution()
    model_files = ['../try_load_model/autoencoder_model_round0.json',
                   '../try_load_model/autoencoder_model_round1.json',
                   '../try_load_model/autoencoder_model_round2.json',
                   '../try_load_model/autoencoder_model_round3.json',
                   '../try_load_model/autoencoder_model_round4.json',
                   '../try_load_model/autoencoder_model_round5.json',
                   '../try_load_model/train_all/autoencoder_model_round6.json']
    weight_files = ['../try_load_model/best_model_round0.h5',
                    '../try_load_model/best_model_round1.h5',
                    '../try_load_model/best_model_round2.h5',
                    '../try_load_model/best_model_round3.h5',
                    '../try_load_model/best_model_round4.h5',
                    '../try_load_model/best_model_round5.h5',
                    '../try_load_model/train_all/best_model_round6.h5']
    train_files = ['../try_load_model/round0_train.dat',
                   '../try_load_model/round1_train.dat',
                   '../try_load_model/round2_train.dat',
                   '../try_load_model/round3_train.dat',
                   '../try_load_model/round4_train.dat',
                   '../try_load_model/round5_train.dat',
                   '../try_load_model/train_all/all_train.dat']
    traj_files = ['output/round_all_weight.dat']
    models = [load_model(x, y) for x, y in zip(model_files, weight_files)]
    train_sets = [pd.read_csv(input_file, delimiter='\s+', comment='#', header=None) for input_file in train_files]
    datasets = [*train_sets, *[read_colvars_traj(x) for x in traj_files]]
    print(datasets[-1])
    models_fve = list()
    for i, dataset in enumerate(datasets):
        dataset_fve = list()
        for j, model in enumerate(models):
            fve = FVE(model, dataset)
            print(f'FVE of model {j} with respect to trajectory {i}: {fve:12.7f}')
            dataset_fve.append(FVE(model, dataset))
        models_fve.append(dataset_fve.copy())
    models_fve = np.array(models_fve)
    plot_fve(models_fve=models_fve, output_filename='EV.png')