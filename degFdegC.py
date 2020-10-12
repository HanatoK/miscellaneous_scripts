#!/usr/bin/env python3
import argparse

def degFTodegC(temp):
    return (temp - 32.0) * 5.0 / 9.0


def degCTodegF(temp):
    return temp * 9.0 / 5.0 + 32.0


def main():
    # parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--FtoC', nargs='+', type=float, help='convert fahrenheit to celsius', default=argparse.SUPPRESS)
    parser.add_argument('--CtoF', nargs='+', type=float, help='convert celsius to fahrenheit', default=argparse.SUPPRESS)
    args = parser.parse_args()

    if 'FtoC' in args:
        for temperature in args.FtoC:
            celsius = degFTodegC(temperature)
            print(f'Convert {temperature:.3f} fahrenheit to celsius: {celsius:.3f}')
    if 'CtoF' in args:
        for temperature in args.CtoF:
            fahrenheit = degCTodegF(temperature)
            print(f'Convert {temperature:.3f} celsius to fahrenheit: {fahrenheit:.3f}')

if __name__ == "__main__":
    # execute only if run as a script
    main()
