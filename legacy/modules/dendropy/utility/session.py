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
Wraps external process as a session, i.e., allow for non-blocking
read/writes to stdout/stderr/stdin.
"""

import sys
import subprocess
import threading

try:
    from Queue import Queue, Empty
except ImportError:
    from queue import Queue, Empty  # python 3.x

ON_POSIX = 'posix' in sys.builtin_module_names

class SessionReader(object):

    def __init__(self, file_handle):
        self.queue = Queue()
        self.stream = file_handle
        self.thread = threading.Thread(
                target=self.enqueue_stream,
                )
        self.thread.daemon = True
        self.thread.start()

    def enqueue_stream(self):
        # for line in self.stream.readline():
        for line in iter(self.stream.readline, b''):
            self.queue.put(line)
        self.stream.close()

    def read(self):
        # read line without blocking
        try:
            line = self.queue.get_nowait()
            # line = self.queue.get(timeout=0.1)
        except Empty:
            return None
        else:
            return line # got line

class Session(object):

    def __init__(self, join_err_to_out=False):
        self.process = None
        self.stdin = None
        self._stdout_reader = None
        self._stderr_reader = None
        self.queue = None
        self.thread = None
        self.join_err_to_out = join_err_to_out

    def start(self, command):
        if self.join_err_to_out:
            stderr = subprocess.STDOUT
        else:
            stderr = subprocess.PIPE
        self.process = subprocess.Popen(command,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=stderr,
                bufsize=1,
                close_fds=ON_POSIX)
        self._stdout_reader = SessionReader(self.process.stdout)
        if not self.join_err_to_out:
            self._stderr_reader = SessionReader(self.process.stderr)

    def _stdin_write(self, command):
        self.process.stdin.write(command)
        self.process.stdin.flush()


