#!/usr/bin/python

"""
detangle.py - Copyright (c) 2012, Howard C. Shaw III
Licensed under the GNU GPL v3

detangle.py [filename] ...

Pass any number of filenames, detangle will extract all trees, and optimize them all simultaneously,
minimizing on all combinations of trees.

Results are output in result1.dat and result2.dat - the program alternates the filename so if you use
ctrl-c to break the execution, one of the two files will have valid data, even if you break during a write.

Results are not the same with every run - the algorithm is stochastic, and may find different local
minima, so multiple runs are recommended.
"""

from sys import argv
import sys
import fileinput
from collections import deque
import itertools
import random
import types

""" Tweak these values to change behavior:

Starting Intensity = the number of twists to try per tree
Number of Iterations = how long after our last improvement do we wait before cooling the process
Max Count = stop after this many iterations - you can use ctrl-c to stop prematurely
Max Iterations = reduce this if you want it to stop when it reaches Intensity 1 and no improvement
is seen for this many iterations
Intensity Reduction = what percent to reduce the intensity to each step
Skip First Tree = 0 for reordering all trees, 1 to leave the first tree fixed

Note that because the iterations are done for each of the trees, the max_iterations may stop
the process before the max_count implies it should, because it is actually counting num_trees times
as fast.
"""

g_starting_intensity = 50
g_number_of_iterations_before_reducing_intensity = 50
g_max_count = 2500
g_max_iterations_without_improvement = 2500
g_intensity_reduction = 0.99
g_skip_first_tree = False
g_output_filename = "result.dat"
g_verbose=False

class tree:
    """ tree: this class encapsulates the individual trees, and anchors the root
         The structure of the tree is never changed once it is created. Instead,
         twists are stored separately and applied on the fly in the output functions.
    """
    def __init__(self, line=None):
        self.root = None
        if not line is None:
            self.name = line[5:line.find(' ',6)]
            self.has_brlens=False
            self.parse(line[line.find('(',6)+1:].strip(';'))
            self.twist_apply_list = self.non_leaves()
            
        #print self.twist_apply_list

    def parse(self, line):
        curstr = line
        self.root = node();
        cur = self.root;
        stack = [self.root];
        if g_verbose:
	        print line
        while len(curstr) > 0:
            if curstr[0] == '(':
                temp = node()
                cur.add(temp)
                stack.append(cur)
                cur = temp
                curstr = curstr[1:]
            elif curstr[0] == ',':
                curstr = curstr[1:]
            elif curstr[0] == ')':
                curstr = curstr[1:]
                if curstr[0]==":":
                    self.has_brlens=True
                    comma = curstr.find(',')
                    paren = curstr.find(')')
                    if comma==-1:
                    	comma=paren
                    if paren==-1:
                    	paren=comma
                    brlen=curstr[1:min(comma,paren)]
                    newstart=min(comma,paren)
                    cur.brlen=float(brlen)
                cur = stack.pop()
                curstr = curstr[newstart:]
#            elif curstr[0] == ':':
#                cur = stack.pop()
#                curstr = curstr[1:]
            else:
                comma = curstr.find(',')
                paren = curstr.find(')')
                colon = curstr.find(':')
                foundcolon=False
                if comma == -1:
                    comma = paren
                if paren == -1:
                    paren = comma
                if colon==-1:
                    colon=comma
                if colon < min(comma,paren):
                    foundcolon=True
                    self.has_brlens=True
                
                if min(comma,paren,colon) > -1:
                    name = curstr[0:min(comma,paren,colon)]
                    if foundcolon:
                        brlen=float(curstr[colon+1:min(comma,paren)])
                    else:
                    	brlen=float(0.0)
                    	
                    if colon!=0:
	                    cur.add(node(name, brlen))
                    curstr = curstr[min(comma,paren):]
                else:
                    curstr = ''
   
    def init_from_phylo(self, phylo):
        self.name = phylo.name
        self.root = node()
        self.init_from_phylo_clade(self.root, phylo.clade)
        self.twist_apply_list = self.non_leaves()

    def init_from_phylo_clade(self, cur, clade):
        if not clade.name is None:
            cur.name = clade.name
        for i in clade.clades:
            n = node()
            self.init_from_phylo_clade(n, i)
            cur.add(n)
            
    def leaves(self):
        d = deque()
        self.root.leaves(d)
        return d

    def non_leaves(self):
        d = deque()
        self.root.non_leaves(d)
        return d        

    def apply_twists(self, twists):
        for i in range(0,min(len(self.twist_apply_list),len(twists))):
            self.twist_apply_list[i].set_twist(twists[i])

    def get_twists(self):
        return [x.twist for x in self.twist_apply_list]

    def print_tree(self):
        print self.name
        self.root.output(1)

    def write(self, f=None):
    	if not f==None:
            f.write("tree ")
            f.write(self.name)
            f.write(" = [&U] ")
            f.write(self.writable())
            f.write(";\n")
        else:
            return self.writable()+";"
        pass

    def writable(self):
        return self.root.writable(self.has_brlens)
    

    def max_depth(self):
        current = self.root
        depth = 0
        max_depth = 0
        stack = [(self.root,0)]
        while len(stack) > 0:
            (current, depth) = stack.pop()
            if current.has_children():
                for i in current.children:
                    max_depth = max(max_depth, depth+1)
                    stack.append((i, depth+1))
        return max_depth
                
    
class node:
    """ node: this class represents all the nodes of a tree, both the branching nodes,
         and the leaf nodes. It will either contain children, or a name.
         It remembers a twist value, and applies it as a rotation to the list of children
         when returning it.
    """
    def __init__(self, name=None, brlen=0.0):
        self.children = []
        self.twist = 0
        self.name = name
        self.brlen=brlen

    def set_twist(self, n):
        self.twist = n

    def add(self, n):
        self.children.append(n)

    def get_children(self):
        d = deque(self.children)
        d.rotate(self.twist)
        return d

    def has_children(self):
        return (len(self.children) > 0)

    def leaves(self, d):
        if len(self.children) == 0:
            d.append(self.name)
        else:
            c = self.get_children()
            for n in c:
                n.leaves(d)

    def non_leaves(self, d):
        if len(self.children) > 0:
            d.append(self)
            c = self.get_children()
            for n in c:
                n.non_leaves(d)

    def output(self, depth):
        if len(self.children) == 0:
            print depth * ' ', '|', self.name
        else:
            d = self.get_children()
            for n in d:
                n.output(depth+1)               

    def writable(self, has_brlens):
        if len(self.children) == 0:
            if has_brlens:
            	return self.name+":"+str(self.brlen)
            else:
                return self.name
        else:
            ret = "("
            d = self.get_children()
            for n in d:
                ret = ret + n.writable(has_brlens) + ","
            ret = ret.rstrip(",") + ")"
            if has_brlens:
            	ret=ret+":"+str(self.brlen)
            return ret
            

def tangle_count_all(trees):
    """ This function applies a tangle counting function to every combination of trees
    """
    l=[]
    for tr in trees.itervalues():
    	leaves=[]
    	for leaf in tr.tree_twisted_leaves():
    		leaves.append(leaf)
    	l.append(leaves)
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
        try:
            if t[a[i]] < t[a[i-1]]:
                count += 1
        except KeyError, e:
            pass
    return count

def flatness_count_all(trees):
    """ This function computes a penalty based on the angle of the lines"""
    l=[]
    for tr in trees.itervalues():
    	leaves=[]
    	for leaf in tr.tree_twisted_leaves():
    		leaves.append(leaf)
    	l.append(leaves)
#    l = [tr.leaves for tr in trees.itervalues()]
    count = 0
    for (a,b) in itertools.combinations(l,2):
        count = count + flatness_count(a,b)
        
    return count

def flatness_count(a, b):
    """ This function computes an angle penalty """
    count = 0
    
    t = dict((b[i],i) for i in range(0,len(b)))
    
    for i in range(0,min(len(a),len(b))):
        try:
            count += abs(i-t[a[i]])
        except KeyError, e:
            pass
              
    return count

def alpha_count_all(trees):
    """ This function computes a penalty for mis-alphabetized trees
    """
    l=[]
    for tr in trees.itervalues():
    	leaves=[]
    	for leaf in tr.tree_twisted_leaves():
    		leaves.append(leaf)
    	l.append(leaves)
    count = 0
    for a in l:
        for i  in range(1,len(a)-1):
            if a[i] < a[i-1]:
                count += 1
    return count

def write(filename, tree_list):
    with open(filename, 'w') as f:
        f.write("#NEXUS \n\n\n")
        f.write("Begin trees;\n")
        for tr in tree_list:
            tr.write(f)
            f.write("\n\n")
        f.write("end;\n")

def minimize_this(trees):
    """ This function needs to return a value to be minimized.
    tangle_count_all counts the crossings
    alpha_count_all counts the alphabetic ordering failures
    you can use either or both, and you can multiply the returned
    values to adjust the importance of alphabetizing vs. tangling
    or you can add your own measure to be minimized. """
    #return tangle_count_all() + (alpha_count_all()*0.5)
    #print flatness_count_all(trees) + tangle_count_all(trees)  + (alpha_count_all(trees)*0.5)
    return flatness_count_all(trees) + tangle_count_all(trees)  + (alpha_count_all(trees)*0.5)




def add_twist_functions(tree):

    
    tree.twist_apply_list=[]
    for node in tree.preorder_internal_node_iter():
        
        node.twist=0
        tree.twist_apply_list.append(node)
        
        def set_twist(self, n):
        	self.twist = n
        node.set_twist = types.MethodType( set_twist, node )
        
        def get_twisted_children(self):
            child_list=[]
            for child in self.child_nodes():
            	child_list.append(child)
            d = deque(child_list)
            d.rotate(self.twist)
            return d
        
        node.get_twisted_children = types.MethodType( get_twisted_children, node )
        
        
        def rotate_twists(self):
            for i in range(0,self.twist):
                self._child_nodes.reverse()
        
        node.rotate_twists = types.MethodType( rotate_twists, node )
    
    for x, node in enumerate(tree.postorder_node_iter()):
        
        def twisted_leaves(self, d):   
            if self.is_leaf():
                d.append(self.taxon.label)
            else:
                c = self.get_twisted_children()
                for n in c:
                    n.twisted_leaves(d)
        
        node.twisted_leaves = types.MethodType( twisted_leaves, node )
        

	
    def get_twists(self):
        return [x.twist for x in self.twist_apply_list]
	
    tree.get_twists = types.MethodType( get_twists, tree )
	
    def apply_twists(self, twists):
        for i in range(0,min(len(self.twist_apply_list),len(twists))):
            self.twist_apply_list[i].set_twist(twists[i])
    
    tree.apply_twists = types.MethodType( apply_twists, tree )
    
    def tree_twisted_leaves(self):
        d = deque()
        self.seed_node.twisted_leaves(d)
        return d
        
    tree.tree_twisted_leaves = types.MethodType( tree_twisted_leaves, tree )
    
    def rotate_nodes(self):
        for node in self.postorder_internal_node_iter():
        	node.rotate_twists()
    
    tree.rotate_nodes = types.MethodType( rotate_nodes, tree )
    
    return tree



tree_list = []
g_starting_intensity = 50
g_number_of_iterations_before_reducing_intensity = 50
g_max_count = 5000
g_max_iterations_without_improvement = 5000
g_intensity_reduction = 0.99
g_skip_first_tree = False

def process_trees(tree_list, starting_intensity=g_starting_intensity,
    number_of_iterations_before_reducing_intensity=g_number_of_iterations_before_reducing_intensity,
    max_count = g_max_count,
    max_iterations_without_improvement = g_max_iterations_without_improvement,
    intensity_reduction = g_intensity_reduction,
    skip_first_tree = g_skip_first_tree,
    output_filename = g_output_filename, verbose = g_verbose):
    """Calculate an initial minimization function value,
    then iteratively take each tree in turn,
    apply _intensity_ random twists to it, and compare the
    overall result with *all* trees for the minimization function.
    Slowly reduce the intensity over time as continued operation
    at a given intensity level ceases to produce improvement.
    """
    first_tree = None
    trees = {}
    twists = {}
    
    for x, tr in enumerate(tree_list):
        if first_tree == None:
            first_tree = tr.label
        if verbose:
            print tr
            
        trees[tr.label] = tr
        tr=add_twist_functions(tr)
        tree_list[x]=tr
        twists[tr.label] = tr.get_twists()
   
    
    if verbose:
	    print tangle_count_all(trees)
    #write(output_filename,tree_list)

    best = minimize_this(trees)
    count = 1
    intensity = starting_intensity
    last_success = 0
    flip = 1
    while intensity > 0 and count < max_count:
    	if verbose:
	        print "Iteration " + str(count) + ", Intensity " + str(intensity) + ", Optimize " + str(best) + ", Tangle Count " + str(tangle_count_all(trees))
        for i in range(0,len(trees)):
            if not skip_first_tree or trees[trees.keys()[i]].label <> first_tree:
                t = list(twists[twists.keys()[i]])
                t2 = list(t)
                for j in range(0,intensity):
                    t[random.randint(0,len(t)-1)] += 1
                trees[trees.keys()[i]].apply_twists(t)
                cur = minimize_this(trees)
                if cur < best:
                    """ If we succeeded in finding a better result, preserve it """
                    last_success = 0
                    twists[twists.keys()[i]] = t
                    best = cur
                    #write("result" + str(flip) + ".dat",tree_list)
                    flip = 3 - flip 
                else:
                    """ Our new result is no better, keep the old tree """
                    trees[trees.keys()[i]].apply_twists(t2)
                    last_success += 1
                if last_success > number_of_iterations_before_reducing_intensity and intensity > 1:
                    intensity = int(intensity * intensity_reduction)
                    last_success = 0
                if last_success > max_iterations_without_improvement and intensity == 1:
                    intensity = 0
        count += 1
    

    for i in range(0,len(tree_list)):
        if not skip_first_tree or i > 0:
#            prelist=[]
#            for node in tree_list[i].leaf_iter():
#            	prelist.append(node.taxon.label)
#            print prelist
            tree_list[i].rotate_nodes()
#            prelist=[]
#            for node in tree_list[i].leaf_iter():
#            	prelist.append(node.taxon.label)
#            print prelist
        
    return (tree_list[0], tree_list[1])
    
    
#    write(output_filename,tree_list)

if __name__=='__main__':
    """
    Loop over all files, reading in all available trees.
    Print the results (for testing/validation purposes).
    """
    for line in fileinput.input():
        #line = line.trim()
        if line[0:4] == 'tree':
            tr = tree(line)
	    if g_verbose:
	            print tr
            tree_list.append(tr)
            #trees[tr.name] = tr
            #twists[tr.name] = tr.get_twists()
            #print trees
    process_trees(tree_list)
                    
                
