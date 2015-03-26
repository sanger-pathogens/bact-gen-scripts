import os, sys

readnames={}
for line in open(sys.argv[1], "rU"):
	words=line.split()
	readnames[words[1]]=int(words[8])


lines=open(sys.argv[2], 'rU').read().split("\n")

output=open(sys.argv[3], "w")

for x in range(0, len(lines), 4):
	readname=lines[x].strip()[1:]
	if len(lines[x].strip())==0:
		continue
	breakpoint=readnames[readname]
	startbit=[]
	endbit=[]
	startbit.append(lines[x]+".1")
	endbit.append(lines[x]+".2")
	startbit.append(lines[x+1][:breakpoint])
	endbit.append(lines[x+1][breakpoint:])
	startbit.append(lines[x+2])
	endbit.append(lines[x+2])
	startbit.append(lines[x+3][:breakpoint])
	endbit.append(lines[x+3][breakpoint:])
	if len(startbit[1])>100:
		print >> output, '\n'.join(startbit)
	if len(endbit[1])>100:
		print >> output, '\n'.join(endbit)

output.close()
