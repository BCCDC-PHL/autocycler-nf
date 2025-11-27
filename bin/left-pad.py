#!/usr/bin/env python

import argparse


def main(args):
    padded_nums = []
    for n in range(1, args.num + 1):
        padded_num = f'{n:02}'
        padded_nums.append(padded_num)

    print(' '.join(padded_nums))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--num', type=int)
    args = parser.parse_args()
    main(args)
