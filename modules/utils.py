import os
import sys
import subprocess
import shlex

class Error (Exception): pass

def open_file_read(filename):
    if filename == '-':
        f = sys.stdin
    elif filename.endswith('.gz') or filename.endswith('.bgz'):
        if not os.path.exists(filename):
            raise Error("Error opening for reading gzipped file '" + filename + "'")

        try:
            f = os.popen('zcat ' + filename)
        except IOError:
            raise Error("Error opening for reading gzipped file '" + filename + "'")
    elif filename.endswith('.bam'):
        try:
            f = os.popen('samtools view -h ' + filename)
        except IOError:
            raise Error("Error opening for reading BAM file '" + filename + "'")
    else:
        try:
            f = open(filename)
        except IOError:
            raise Error("Error opening for reading file '" + filename + "'")

    return f

def open_file_write(filename):
    if filename == '-':
        f = sys.stdout
    elif filename.endswith('.bgz'):
        if not os.path.exists(os.path.abspath(os.path.dirname(filename))):
            raise Error("Error opening for writing bgzipped file '" + filename + "'")
        try:
            f = os.popen('bgzip -c > ' + filename, 'w')
        except IOError:
            raise Error("Error opening for writing bgzipped file '" + filename + "'")
    elif filename.endswith('.gz'):
        if not os.path.exists(os.path.abspath(os.path.dirname(filename))):
            raise Error("Error opening for writing gzipped file '" + filename + "'")

        try:
            f = os.popen('gzip -9 -c > ' + filename, 'w')
        except IOError:
            raise Error("Error opening for writing gzipped file '" + filename + "'")
    else:
        try:
            f = open(filename, 'w')
        except IOError:
            raise Error("Error opening for writing file '" + filename + "'")

    return f

def close(filehandle):
    if filehandle not in [sys.stdout, sys.stderr]:
        filehandle.close()

def file_column_combine(action_list, file_list, outfile):
    SUM = 0
    MAX = 1
    MIN = 2
    SAME = 3

    convert_actions = {'sum': SUM,
                       'max': MAX,
                       'min': MIN,
                       'same': SAME}
    try:
        actions = [convert_actions[x] for x in action_list]
    except KeyError:
        raise Error('Not recognised in file_column_combine:\n' + str(action_list))

    file_handles = []
    fout = open_file_write(outfile)
    number_indices = [i for i in range(len(action_list)) if action_list[i] in ['min', 'max', 'sum']]

    for fname in file_list:
        file_handles.append(open_file_read(fname))

    for line in file_handles[0]:
        lines = [line.rstrip().split()]
        for fh in file_handles[1:]:
            lines.append(fh.readline().rstrip().split())
            if len(lines[0]) != len(lines[-1]):
                raise Error('Mismatch in number of columns\n' + lines)


        if line.startswith('#'):
            print >> fout, line.rstrip()
        else:
            out = []
            for i in range(len(actions)):
                if actions[i] == SAME:
                    out.append(lines[0][i])
                elif actions[i] == MAX:
                    out.append(max([float(a[i]) for a in lines]))
                elif actions[i] == MIN:
                    out.append(min([float(a[i]) for a in lines]))
                elif actions[i] == SUM:
                    out.append(sum([float(a[i]) for a in lines]))

            for i in number_indices:
                if out[i].is_integer():
                    out[i] = int(out[i])
            print >> fout, '\t'.join([str(x) for x in out])

    for fh in file_handles:
        close(fh)

    close(fout)


def file_transpose(f_in, f_out, sep_in=None, sep_out='\t'):
    rows = []
    f = open_file_read(f_in)
    for line in f:
        rows.append(line.rstrip().split(sep_in))
    close(f)

    columns_out = max([len(x) for x in rows])

    for r in rows:
        r += ['.'] * (columns_out - len(r))

    f = open_file_write(f_out)
    for i in range(columns_out):
        print >> f, sep_out.join([str(rows[x][i]) for x in range(len(rows))])

    close(f)

def syscall(cmd):
    retcode = subprocess.call(cmd, shell=True)

    if retcode != 0:
        raise Error("Error in system call. Command was:\n" + cmd)

def syscall_get_stdout(cmd):
    try:
        out = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE).communicate()[0].decode('utf-8').rstrip()
        return out.split('\n')
    except:
        raise Error('Error in system call. I tried to run:\n' + str(cmd))


