#!/usr/bin/python

"""
dendropy_detangler.py - Copyright (c) 2012, Howard C. Shaw III
Licensed under the GNU GPL v3

dendropy_detangler.py [filename] ...

Pass any number of filenames, detangler will extract all trees using dendropy,
and optimize them all simultaneously, minimizing on all combinations of trees.

WARNING: dendropy appears to strip the names under some conditions,
which can break the detangle algorithm. Specifically, on godef.tre, with names of the
form PAUP_1, it returned the name PAUP for all trees.

Also, dendropy appears to choke on Nexus files containing UTF8 Byte-Order-Marks.
Use this
tail --bytes=+4 test.dat > newfile.dat
to strip them.
"""

from sys import argv
import dendropy
from StringIO import StringIO
from detangle import tree, node, process_trees
import argparse

if __name__=='__main__':
    parser = argparse.ArgumentParser(description = 'Minimize tangling across multiple trees.')
    #parser.add_argument('-t, --output-type', default='nexus', choices=['newick', 'nexus', 'phyloxml'])
    parser.add_argument('-o, --output-filename', dest='output_filename', nargs='?', default='result.dat')
    parser.add_argument('infiles', nargs='+')
    args = parser.parse_args()
    print args
    dpy_trees = {}
    dpy_tree_list = []
    for filename in args.infiles:
        """ Work out what format the file is, then parse the trees out """
        f = open(filename, 'r')
        first = f.readline()
        if first.find('#nexus') > -1 or first.find('#NEXUS') > -1:
            tree_type='nexus'
        elif first.find('<') > -1:
            tree_type='phyloxml'
        else:
            tree_type='newick'
        f.close()

        if tree_type == 'newick':
            # only one tree per file, so get a tree
            tree = dendropy.Tree.get_from_path(filename, tree_type)
            dpy_trees[tree.label] = tree
            dpy_tree_list.append(tree)
        else:
            # the others can have more than one per,
            # so get a TreeList
            t = dendropy.TreeList.get_from_path(filename, tree_type)
            if not t is None:
                for i in t:
                    dpy_trees[i.label] = i
                    dpy_tree_list.append(i)
    trees = {}
    tree_list = [ ]
    """ At this point, dpy_trees is full, now we need to convert these dendropy trees to the detangle
    trees which are optimized for rotations. Luckily, dendropy trees have a utility function that
    provides us with exactly the simple Newick string detangle.tree likes."""
    for i in dpy_tree_list:
        tr = tree('tree ' + i.label + ' = [&U] ' + i.as_newick_string())
        trees[tr.name]=tr
        tree_list.append(tr)
    process_trees(tree_list, output_filename = args.output_filename)
