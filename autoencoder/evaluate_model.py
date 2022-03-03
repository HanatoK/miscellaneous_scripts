#!/usr/bin/env python3
from tensorflow.keras.models import model_from_json
import tensorflow as tf
import pandas as pd
import argparse


def renormalize_weight(weight):
    renorm_factor = len(weight) / weight.sum()
    return weight * renorm_factor


def compile_model(model):
    model.compile(optimizer=tf.keras.optimizers.SGD(momentum=0.1, nesterov=True), loss='mse', metrics=['mae'])
    return model


def load_model(model_file, weight_file, compile=False):
    with open(model_file, 'r') as json_file:
        model = model_from_json(json_file.read())
        model.load_weights(weight_file, by_name=True)
        if compile:
            model = compile_model(model)
        return model


def evaluate_model(model, traj_filename, num_variables=12):
    data = pd.read_csv(traj_filename, delimiter='\s+', header=None, comment='#')
    w = None
    x = data[list(range(0, num_variables))]
    if data.shape[1] == num_variables+1:
        w = renormalize_weight(data[num_variables])
    results = model.evaluate(x, x, sample_weight=w, batch_size=256)
    return results


def predict_data(model, traj_filename, num_variables=12):
    data = pd.read_csv(traj_filename, delimiter='\s+', header=None, comment='#')
    # w = None
    x = data[list(range(0, num_variables))]
    # if data.shape[1] == num_variables+1:
    #     w = renormalize_weight(data[num_variables])
    results = model.predict(x, batch_size=data.shape[0])
    return results


def is_autoencoder(model):
    firsttime = True
    shape_first_layer = 0
    for layer in model.layers:
        if firsttime:
            shape_first_layer = layer.get_config()['batch_input_shape']
            firsttime = False
        shape_last_layer = layer.output_shape
    if shape_first_layer == shape_last_layer:
        return True
    else:
        # print('The shape of the first layer:')
        # print(shape_first_layer)
        # print('The shape of the last layer:')
        # print(shape_last_layer)
        return False


def get_input_size(model):
    shape_first_layer = model.layers[0].get_config()['batch_input_shape']
    return shape_first_layer[1]


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('model_file', type=str, help='Encoder model file (JSON format)')
    parser.add_argument('weights_file', type=str, help='Saved weights and biases (HDF5 format)')
    parser.add_argument('dataset', type=str, help='Dataset to predict')
    parser.add_argument('output', type=str, help='Output of prediction')
    args = parser.parse_args()
    # load model
    model = load_model(model_file=args.model_file, weight_file=args.weights_file)
    model.summary()
    # predict data
    results = predict_data(model=model, traj_filename=args.dataset, num_variables=get_input_size(model))
    # write to file
    with open(args.output, 'w') as f_output:
        df = pd.DataFrame(results)
        df.to_string(f_output, float_format='%15.10f', index=False, header=False)
        f_output.write('\n')
    # if this is an autoencoder model, evaluate the mse
    if is_autoencoder(model):
        print('This is an AE model.')
        model = compile_model(model)
        print(evaluate_model(model=model, traj_filename=args.dataset, num_variables=get_input_size(model)))
