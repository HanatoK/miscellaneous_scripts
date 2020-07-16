#!/usr/bin/env python3
import numpy as np

def readFepLog(log_filename):
    forward_lambda = []
    forward_delta_A = []
    forward_delta_epsilon = []
    backward_lambda = []
    backward_delta_A = []
    backward_delta_epsilon = []
    with open(log_filename) as fp:
        print(f'Reading {log_filename}')
        for line in fp:
            if line.strip():
                line_data = line.split()
                if line_data[0] == 'forward:':
                    forward_lambda.append(float(line_data[1]))
                    forward_delta_A.append(float(line_data[3]))
                    forward_delta_epsilon.append(float(line_data[4]))
                if line_data[0] == 'backward:':
                    backward_lambda.append(float(line_data[1]))
                    backward_delta_A.append(float(line_data[3]))
                    backward_delta_epsilon.append(float(line_data[4]))
    forward_data = [np.array(forward_lambda), np.array(forward_delta_A), np.array(forward_delta_epsilon)]
    backward_data = [np.array(backward_lambda), np.array(backward_delta_A), np.array(backward_delta_epsilon)]
    return [forward_data, backward_data]

def statAvgVar(arg1, *argv):
    data_sum = arg1
    data_sum2 = arg1 * arg1
    count = 1
    for arg in argv:
        data_sum += arg
        data_sum2 += arg * arg
        count += 1
    average = data_sum / float(count)
    average2 = data_sum2 / float(count)
    variance = average2 - average * average
    return [average, variance]

# arg1: output
# argv: logs from parsefep
def analyzeFepLog(arg1, *argv):
    forward_delta_A_array = []
    backward_delta_A_array = []
    for logfile in argv:
        forward_data, backward_data = readFepLog(logfile)
        forward_delta_A_array.append(forward_data[1])
        backward_delta_A_array.append(backward_data[1])
        forward_lambda = forward_data[0]
        backward_lambda = backward_data[0]
    forward_avg, forward_var = statAvgVar(*forward_delta_A_array)
    backward_avg, backward_var = statAvgVar(*backward_delta_A_array)
    forward_output = arg1 + '.forw.dat'
    backward_output = arg1 + '.back.dat'
    np.savetxt(forward_output, np.c_[forward_lambda, forward_avg, np.sqrt(forward_var)], fmt = '%12.7f')
    np.savetxt(backward_output, np.c_[backward_lambda, backward_avg, np.sqrt(backward_var)], fmt = '%12.7f')
    return

input_file_list = []
for i in range(1,6):
    log_filename = '/home/yjcoshc/HDD/md/FEP-test/01.ethane-ethane/' + 'pme_gpu_' + str(i) + '/ParseFEP.log'
    input_file_list.append(log_filename)
analyzeFepLog('pme_gpu', *input_file_list)
