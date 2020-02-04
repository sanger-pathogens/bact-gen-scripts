import socket
import utils
import subprocess

class Error (Exception): pass

class Bsub:
    """Class for submitting bsub jobs on the farm"""
    def __init__(self, out, error, name, queue, mem, cmd,
                 start=0, end=0,
                 more="",
                 depend=None,
                 threads=0,
                 tmp_space=0,
                 array_list=None,
                 ended=None,
                 max_array_size=100):

        self.q = queue
        self.m = round(mem, 3)
        self.extra = more
        self.c = cmd
        self.threads = threads
        self.tmp_space = tmp_space

        if self.threads > 1:
            self.extra += ' -n' + str(threads)

        # is this a job array?
        if 0 < start:
            self.o = out + '.%I'
            self.e = error + '.%I'
            self.n = '"' + name + '[' + str(start) + '-' + str(end) + ']%' + str(max_array_size) + '"'

            # translate every appearance of 'INDEX' in command to \$LSB_JOBINDEX
            self.c = self.c.replace('INDEX', '\$LSB_JOBINDEX')
        elif array_list is not None:
            self.o = out + '.%I'
            self.e = error + '.%I'
            self.n = '"' + name + '[' + ','.join(array_list) + ']%' + str(max_array_size) + '"'
            self.c = self.c.replace('INDEX', '\$LSB_JOBINDEX')
        else:
            self.o = out
            self.e = error
            self.n = name

        self.run_when_done = []
        self.run_when_ended = []
        self.add_dependency(depend)
        self.add_dependency(ended, ended=True)


    def __str__(self):
        # work out the memory in megs/kb
        mb = int(self.m * 1000)
        kb = int(str(mb) + '000')
        tmp_space_select = ''
        tmp_space_rusage = ''

        if len(self.run_when_done) + len(self.run_when_ended) > 0:
            dependencies =  ' && '.join(['done(' + x + ')' for x in self.run_when_done] + ['ended(' + x + ')' for x in self.run_when_ended])
            dependencies = "-w'" + dependencies + "'"
        else:
            dependencies = ''

        syscmd = "bsub -E 'test -e /nfs/users/nfs_s/sh16/ '" + self.extra + " " + dependencies + " -q " + self.q + ' -R "'
        if self.threads > 1:
            syscmd += 'span[hosts=1] '

        syscmd += 'select[type==X86_64 && mem > ' + str(mb)
        if self.tmp_space:
            tmp_mb = int(self.tmp_space * 1000)
            syscmd += ' && tmp>' + str(tmp_mb)

        syscmd += '] rusage[mem=' + str(mb)


        if self.tmp_space:
            syscmd += ', tmp=' + str(tmp_mb)

        # New farm specifies all memory in MB. Other farms use kb for the -M option
        cmd = 'lsadmin showconf lim ' + socket.gethostname() + r''' | grep LSF_UNIT_FOR_LIMITS || echo 'LSF_UNIT_FOR_LIMITS = KB' '''
        lsf_units = subprocess.check_output(cmd, shell=True).decode('utf-8').rstrip().split()[2]

        if lsf_units == 'MB':
            syscmd += ']" -M' + str(mb)
        elif lsf_units == 'KB':
            syscmd += ']" -M' + str(kb)
        else:
            raise Error('Error getting lsf memory units. cannot continue')

        syscmd += ' -o ' + self.o + ' -e ' + self.e + ' -J ' + self.n + ' ' + self.c
        return syscmd

    # takes a list of jod ids, makes this object depend on them to run.
    # If ended = true, makes this job run when the given job ends, regardless of its success.
    # These are equivalent to -w'done(...)' and -w'ended(...)'
    def add_dependency(self, deps, ended=False):
        if deps is None or len(deps) == 0:
            return

        if type(deps) is not list:
            deps = [deps]

        for d in deps:
            try:
                x = int(d)
                x = d
            except ValueError:
                x = '"' + d + '"'

            if ended:
                self.run_when_ended.append(x)
            else:
                self.run_when_done.append(x)

    # submits the job and returns the job ID if submission successful.  Dies otherwise
    def run(self, verbose=False):
        if verbose:
            print(self)

        try:
            # for ok bsub call, output should look like this:
            # Job <42> is submitted to queue <normal>.
            b_out = utils.syscall_get_stdout(str(self))[0]
        except:
            raise Error('Error in bsub call. I tried to run:\n' + str(self))

        if 'is submitted to' in b_out:
            if verbose:
                print(b_out)

            return b_out.split()[1][1:-1]
        else:
            raise Error('Error in bsub call. I tried to run:\n' + str(self))


class JobStats:
    def __init__(self, bsub_stdout_file):
        self.name = 'test'
