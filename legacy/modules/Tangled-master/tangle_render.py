#!/usr/bin/python

"""
tangle_render.py - Copyright (c) 2012, Howard C. Shaw III
Licensed under the GNU GPL v3

tangle_render.py [filename] ...

Pass any number of filenames, detangle will extract all trees, and render tangles on all
combinations of trees.
"""

line_gap = 5
line_region_width = 260
line_darkness = 0.3

from detangle import tree, node
from sys import argv
import fileinput
from collections import deque
import itertools
import random
import time

import cairo
import rsvg

def tangle_count_all():
    """ This function applies a tangle counting function to every combination of trees
    """
    l = [tr.leaves() for tr in trees.itervalues()]
    count = 0
    for (a,b) in itertools.combinations(l,2):
        count = count + tangle_count(a,b)
    return count
    
def tangle_count(a, b):
    """ This function computes a tangle count by counting the number of times a
         pair of leaves in the left tree are in the opposite order in the right tree
    """
    count = 0
    t = dict((b[i],i) for i in range(0,len(b)))
    
    for i in range(1,min(len(a),len(b))):
        if t[a[i]] > t[a[i-1]]:
            count += 1
    return count

def alpha_count_all():
    """ This function computes a penalty for mis-alphabetized trees
    """
    l = [tr.leaves() for tr in trees.itervalues()]
    count = 0
    for a in l:
        for i  in range(1,len(a)-1):
            if a[i] < a[i-1]:
                count += 1
    return count

def write(filename):
    with open(filename, 'w') as f:
        f.write("#NEXUS \n\n\n")
        f.write("Begin trees;\n")
        for tr in trees.itervalues():
            tr.write(f)
            f.write("\n\n")
        f.write("end;\n")

tree_list = []        
trees = {}
twists = {}
first_tree = None

def draw_tree(ct, tr, x_pos, height):
    ct.set_source_rgb(0, 0, 0)
    widest = 0
    x, y, w, h = ct.text_extents("Tree: " + tr.name)[:4]
    widest = w
    ct.move_to(x_pos, 0)
    ct.show_text("Tree: " + tr.name)
    vert = height
    for l in tr.leaves():
        x, y, w, h = ct.text_extents(l)[:4]
        widest = max(w, widest)
        ct.move_to(x_pos, vert)
        ct.show_text(l)
        vert += height
    return widest

def draw_lines(ct, left, right, x1, x2, height):
    ct.set_source_rgba(0, 0, 0, line_darkness)
    b = right.leaves()
    t = dict((b[i],i) for i in range(0,len(b)))
    i = 1
    for l in left.leaves():
        if t.has_key(l):
            ct.move_to(x1, height * i - (height / 3))
            ct.line_to(x2, height * (t[l] + 1) - (height / 3))
            ct.stroke()
        i += 1

def draw_branches(ct, tr, x_pos, y_pos, height, branch_length, facing_right = True, depth = 0, current = None, max_depth = 0):   
    if current == None:
        current = tr.root
        max_depth = tr.max_depth()
    top = 0.00
    upper = None
    lower = 0.00
    total_height = 0
    middle_height = 0
    #print y_pos
    children = current.get_children()
    for i in children:
        if (i.has_children()):
            (total_height, middle_height) = draw_branches(ct, tr, x_pos, y_pos + top, height, branch_length, facing_right, depth + 1, i, max_depth)
            if facing_right:
                ct.move_to(x_pos + (branch_length * depth), y_pos + top + middle_height)
                ct.line_to(x_pos + (branch_length * (depth + 1)), y_pos + top + middle_height)
            else:
                ct.move_to(x_pos + (branch_length * (max_depth-depth)), y_pos + top + middle_height)
                ct.line_to(x_pos + (branch_length * ((max_depth-depth) - 1)), y_pos + top + middle_height)
            ct.stroke()
            if (upper == None):
                upper = top + middle_height
            lower = max(lower, top + middle_height)
            top += total_height
        else:
            if (upper == None):
                upper = top
            if facing_right:
                ct.move_to(x_pos + (branch_length * depth), y_pos + top)
                ct.line_to(x_pos + (branch_length * max_depth), y_pos + top )
            else:
                ct.move_to(x_pos + (branch_length * (max_depth-depth)), y_pos + top)
                ct.line_to(x_pos + 0, y_pos + top )
            ct.stroke()
            lower = max(lower, top)
            top += height
    if facing_right:
        ct.move_to(x_pos + (branch_length * depth), y_pos + upper)
        ct.line_to(x_pos  + (branch_length * depth), y_pos + lower)
    else:
        ct.move_to(x_pos + (branch_length * (max_depth-depth)), y_pos + upper)
        ct.line_to(x_pos  + (branch_length * (max_depth-depth)), y_pos + lower)
    ct.stroke()
    if current == tr.root:
        return (max_depth+1) * branch_length
    else:
        return (top, upper + ((lower - upper)/2.00))

if __name__=='__main__':
    """
    Loop over all files, reading in all available trees.
    """
    for line in fileinput.input():
        #line = line.trim()
        if line[0:4] == 'tree':
            tr = tree(line)
            if first_tree == None:
                first_tree = tr.name
            trees[tr.name] = tr
            tree_list.append(tr)
            twists[tr.name] = tr.get_twists()
            #print trees

    combos = list(itertools.combinations(tree_list,2))

    f = open('output.svg', 'w')
    surf = cairo.SVGSurface(f, len(combos) * (120 + line_gap+ line_region_width) * 2, len(tr.leaves()) * 40)
    ct = cairo.Context(surf)
    ct.translate(10,16)
    ct.set_source_rgb(0.0, 0.0, 0.0)
    ct.select_font_face("Georgia", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
    ct.set_font_size(12)
    
    x_bearing, y_bearing, width, height = ct.text_extents("Tygp")[:4]
    height = height * 1.25
    left_x = 0
    left_tree = tree_list[0]

    w = draw_tree(ct, left_tree, 0, height)

    while len(combos) > 0:
        valid = list(x for x in combos if left_tree in x)
        #print "Combos " , list(x for x in combos), " Valid " , list((a.name, b.name) for (a,b) in (x for x in valid))
        if len(valid) > 0:
            if valid[0][0] == left_tree:
                right_tree = valid[0][1]
            else:
                right_tree = valid[0][0]
            """ We have a left tree and a right tree - draw the right tree, and the links """
            draw_lines(ct, left_tree, right_tree, left_x + w + line_gap, left_x + w + line_gap + line_region_width, height)
            left_x += w + line_gap + line_region_width + line_gap
            w = draw_tree(ct, right_tree, left_x, height)
            left_tree = right_tree
            combos.remove(valid[0])
        else:
            left_x += w + line_gap
            left_tree = combos[0][0]
            w = draw_tree(ct, left_tree, left_x, height)
    surf.finish()
    #time.sleep(5)
    f.close()
            
    combos = list(itertools.combinations(tree_list,2))
    while len(combos) > 0:
        left_tree = combos[0][0]
        right_tree = combos[0][1]

        f = open(left_tree.name + ':' + right_tree.name + '.svg', 'w')
        surf = cairo.SVGSurface(f, left_tree.max_depth() * 20 + right_tree.max_depth() * 20 + len(combos) * (120 + line_gap+ line_region_width) * 2, len(tr.leaves()) * 40)
        ct = cairo.Context(surf)
        ct.translate(10,16)
        ct.set_source_rgb(0.0, 0.0, 0.0)
        ct.select_font_face("Georgia", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
        ct.set_font_size(12)
        
        x_bearing, y_bearing, width, height = ct.text_extents("Tygp")[:4]
        height = height * 1.25
        left_x = 0
        
        """ We have a left tree and a right tree - draw the right tree, and the links """
        w = draw_branches(ct, left_tree, 0, height - height / 3, height, 20, True)
        left_x = w
        w = draw_tree(ct, left_tree, left_x + line_gap, height)
        draw_lines(ct, left_tree, right_tree, left_x + w + line_gap, left_x + w + line_gap + line_region_width, height)
        left_x += w + line_gap + line_region_width + line_gap
        w = draw_tree(ct, right_tree, left_x, height)
        w = draw_branches(ct, right_tree, left_x + w + line_gap, height - height / 3, height, 20, False)
        combos.remove(combos[0])
        surf.finish()
        #time.sleep(5)
        f.close()
            

                    
                
