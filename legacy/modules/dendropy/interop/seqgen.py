#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010 Jeet Sukumaran and Mark T. Holder.
##  All rights reserved.
##
##  See "LICENSE.txt" for terms and conditions of usage.
##
##  If you use this work or any portion thereof in published work,
##  please cite it as:
##
##     Sukumaran, J. and M. T. Holder. 2010. DendroPy: a Python library
##     for phylogenetic computing. Bioinformatics 26: 1569-1571.
##
##############################################################################

"""
Wrappers for interacting with SeqGen. Originally part of PySeqGen.
"""

import subprocess
from cStringIO import StringIO
import uuid
import tempfile
import socket
import random
import os
import sys

from optparse import OptionGroup
from optparse import OptionParser

import dendropy
from dendropy.utility.messaging import get_logger
_LOG = get_logger("interop.seqgen")

HOSTNAME = socket.gethostname()
PID = os.getpid()

def _get_strongly_unique_tempfile(dir=None):
    return tempfile.NamedTemporaryFile(dir=dir, prefix="dendropy_tempfile-{0}-{1}-{2}".format(HOSTNAME, PID, uuid.uuid4()))

def _get_tempfile(dir=None):
    return tempfile.NamedTemporaryFile(dir=dir)

class SeqGen(object):
    """
    This class wraps all attributes and input needed to make a call to SeqGen.
    """

    class SubstitutionModel(object):
        def __init__(self, idstr):
            self.idstr = idstr
        def __str__(self):
            return self.idstr

    F84 = SubstitutionModel("F84")
    HKY = SubstitutionModel("HKY")
    GTR = SubstitutionModel("GTR")
    JTT = SubstitutionModel("JTT")
    WAG = SubstitutionModel("WAG")
    PAM = SubstitutionModel("PAM")
    BLOSUM = SubstitutionModel("BLOSUM")
    MTREV = SubstitutionModel("MTREV")
    CPREV = SubstitutionModel("CPREV")
    GENERAL = SubstitutionModel("GENERAL")
    MODELS = [F84, HKY, GTR, JTT, WAG, PAM, BLOSUM, MTREV, CPREV, GENERAL]
    MODEL_IDS = [str(m) for m in MODELS]

    def get_model(idstr):
        for model in SeqGen.MODELS:
            if idstr.upper() == model.idstr.upper():
                return model
        return None
    get_model = staticmethod(get_model)

    def __init__(self, strongly_unique_tempfiles=False):
        """
        Sets up all properties, which (generally) map directly to command
        parameters of Seq-Gen.
        """

        # control tempfile generation
        if strongly_unique_tempfiles:
            self.get_tempfile = _get_strongly_unique_tempfile
        else:
            self.get_tempfile = _get_tempfile

        # python object specific attributes
        self.seqgen_path = 'seq-gen'
        self.rng_seed = None
        self._rng = None

        # following are passed to seq-gen in one form or another
        self.char_model = 'HKY'
        self.seq_len = None
        self.num_partitions = None
        self.scale_branch_lens = None
        self.scale_tree_len = None
        self.codon_pos_rates = None
        self.gamma_shape = None
        self.gamma_cats = None
        self.prop_invar = None
        self.state_freqs = None
        self.ti_tv = 0.5 # = kappa of 1.0, i.e. JC
        self.general_rates = None
        self.ancestral_seq = None
        self.output_text_append = None
        self.write_ancestral_seqs = False
        self.write_site_rates = False

    def _get_rng(self):
        if self._rng is None:
            self._rng = random.Random(self.rng_seed)
        return self._rng
    def _set_rng(self, rng):
        self._rng = rng
    rng = property(_get_rng, _set_rng)

    def _get_kappa(self):
        return float(self.ti_tv) / 2
    def _set_kappa(self, kappa):
        self.ti_tv = kappa * 2

    def _compose_arguments(self):
        """
        Composes and returns a list of strings that make up the arguments to a Seq-Gen
        call, based on the attribute values of the object.
        """
        args = []
        args.append(self.seqgen_path)
        if self.char_model:
            args.append("-m%s" % str(self.char_model))
        if self.seq_len:
            args.append("-l%s" % self.seq_len)
        if self.num_partitions:
            args.append("-p%s" % self.num_partitions)
        if self.scale_branch_lens:
            args.append("-s%s" % self.scale_branch_lens)
        if self.scale_tree_len:
            args.append("-d%s" % self.scale_tree_len)
        if self.codon_pos_rates:
            args.append("-c%s" % (",".join(self.codon_pos_rates)))
        if self.gamma_shape:
            args.append("-a%s" % self.gamma_shape)
        if self.gamma_cats:
            args.append("-g%s" % self.gamma_cats)
        if self.prop_invar:
            args.append("-i%s" % self.prop_invar)
        if self.state_freqs:
            if isinstance(self.state_freqs, str):
                args.append("-f%s" % self.state_freqs)
            else:
                args.append("-f%s" % (",".join([str(s) for s in self.state_freqs])))
        if self.ti_tv and (self.char_model in ['HKY', 'F84']):
            args.append("-t%s" % self.ti_tv)
        if self.general_rates:
            if isinstance(self.general_rates, str):
                args.append("-r%s" % self.general_rates)
            else:
                args.append("-r%s" % (",".join([str(r) for r in self.general_rates])))
        if self.ancestral_seq:
            args.append("-k%s" % self.ancestral_seq)
        if self.output_text_append:
            args.append("-x'%s'" % self.output_text_append)
        if self.write_ancestral_seqs:
            args.append("-wa")
        if self.write_site_rates:
            args.append("-wr")

        # following are controlled directly by the wrapper
        # silent running
        args.append("-q")
        # we explicitly pass a random number seed on each call
        args.append("-z%s" % self.rng.randint(0, sys.maxint))
        # force nexus
        args.append("-on")
        # force one dataset at a time
        args.append("-n1")
        return args

    def generate(self, trees, dataset=None, taxon_set=None, **kwargs):
        args=self._compose_arguments()
        tree_inputf = self.get_tempfile()
        trees.write_to_path(tree_inputf.name,
                "newick",
                suppress_rooting=True,
                suppress_internal_node_labels=True)
        tree_inputf.flush()
        args.append(tree_inputf.name)
        #_LOG.debug("seq-gen args: = %s" % " ".join(args))
        run = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = run.communicate()
        if stderr or run.returncode != 0:
            raise RuntimeError("Seq-gen error: %s" % stderr)
        if taxon_set is None:
            taxon_set = trees.taxon_set
        if dataset is None:
            dataset = dendropy.DataSet(taxon_set=taxon_set, **kwargs)
        results = StringIO(stdout)
        #_LOG.debug('stderr = ' + stderr)
        #_LOG.debug('stdout = ' + stdout)

        dataset.read(results, "nexus")
        return dataset



