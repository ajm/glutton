import sys

i = 1

for line in sys.stdin :
    print ">%d\n%s" % (i, line.strip())
    i += 1

