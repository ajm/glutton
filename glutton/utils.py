import tempfile
import os

def tmpfile() :
    fd,name = tempfile.mkstemp(prefix='glutton', dir=os.environ['TMPDIR'])
    os.close(fd)
    return name

