import os
import subprocess

from sys import exit
from os.path import isfile

from glutton.utils import get_log, get_binary_path


class ExternalToolError(Exception) :
    pass

class ExternalTool(object) :
    def __init__(self, location=None) :
        self.binary_location = get_binary_path(self.name) if not location else location
        self.log = get_log()

    @property
    def name(self) :
        return type(self).__name__.lower()

    @property
    def version(self) :
        raise NotImplementedError()

    def _execute(self, parameters, expected_outfiles) :
        returncode = 0
        output = ""

        self.log.debug(' '.join([self.binary_location] + parameters))

        try :
            output = subprocess.check_output(
                                [self.binary_location] + parameters, 
                                stderr=subprocess.STDOUT,
                                close_fds=True
                                )

        except subprocess.CalledProcessError, cpe :
            returncode = cpe.returncode
            output = cpe.output

        # some program misbehave, so be careful to check the expected output
        # and change returncode as necessary
        if returncode == 0 :
            missing = [o for o in expected_outfiles if not isfile(o)]
            if len(missing) != 0 :
                returncode = 256
                output = "the following files were missing : %s" % ' '.join(missing)

        return returncode, output

