#!/usr/bin/env python3
import numpy as np
from tensorflow.keras.models import model_from_json
import argparse


def load_model(model_filename, weights_filename):
    '''load tensorflow model from saved files'''
    with open(model_filename, 'r') as json_file:
        model = model_from_json(json_file.read())
        model.load_weights(weights_filename)
        return model


def dump_dense_weights_biases(model):
    for index, layer in enumerate(model.layers):
        name = layer.name
        data = layer.get_weights()
        if len(data) == 0:
            continue
        if name.startswith('dense'):
            # print(name)
            weights = data[0]
            biases = data[1]
            np.savetxt(f'dense_{index}_weights.txt', weights.transpose())
            np.savetxt(f'dense_{index}_biases.txt', biases)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('model_file', type=str, help='model file of previous round')
    # parser.add_argument('--previous_encoder_model', type=str, help='encoder model of previous round')
    parser.add_argument('weights_file', type=str, help='weight file of previous round')
    args = parser.parse_args()
    model = load_model(args.model_file, args.weights_file)
    dump_dense_weights_biases(model)
