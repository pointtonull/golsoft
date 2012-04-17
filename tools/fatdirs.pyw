#!/usr/bin/env python
#-*- coding: UTF-8 -*-

from math import log
import os
import sys

def main():

    if len(sys.argv) == 2:
        path = sys.argv[1]
    else:
        path = "./"

    weighs = []
    for root, dirs, files in os.walk(path):

        dotfiles = [dotfile for dotfile in files if dotfile.startswith(".")]
        for dotfile in dotfiles:
            files.remove(dotfile)

        dotdirs = [dotdir for dotdir in dirs if dotdir.startswith(".")]
        for dotdir in dotdirs:
            dirs.remove(dotdir)

        if not root.startswith("./."):
            weigh = len(dirs) + log(len(files) + 1)
            weighs.append((weigh, root))

    weighs.sort()
    for weigh, root in weighs[-20:]:
        print("%6.2f %s" % (weigh, root))


if __name__ == "__main__":
    exit(main())
