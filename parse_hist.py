#!/usr/bin/python3
import argparse
import os
from subprocess import call

parser = argparse.ArgumentParser()
parser.add_argument("hist", help = "specify the .hist.grad or .hist.count file")
args = parser.parse_args()
histfilename = os.path.basename(args.hist)
newbasename = os.path.splitext(histfilename)[0] + "_"
extname = os.path.splitext(histfilename)[1]
count = 0
histfile = open(histfilename, "r")
firstline = histfile.readline().rstrip()
for line in histfile:
    count = count + 1
    nline = line.rstrip()
    if nline == firstline:
        print(count)
        break
histfile.close()
extarg = "--additional-suffix=" + extname
digiarg = "--numeric-suffixes=1"
# 如果hist文件分割出来多于9999个，请酌情修改下方的后缀长度
suffixlenarg = "--suffix-length=4"
call(["split", "-l", str(count), args.hist, newbasename, extarg, digiarg, suffixlenarg])