import sys
import os
import re
import datetime

from lib.datatypes import Sequence, IUPAC


class DataFileError(Exception):
    pass

class DataFile(object) :
    def __init__(self, fname, extension) :
        self.pathname = os.path.dirname(fname)
        self.filename = os.path.basename(fname)
        self.extension = extension

        if not self.__exists() :
            raise DataFileError("'%s' does not exist" % self.get_filename())

    @property
    def name(self) :
        return self.get_filename()

    def __exists(self) :
        return os.path.isfile(self.get_filename())

    def get_filename(self) :
        if self.pathname != "" :
            return self.pathname + os.sep + self.filename
        return self.filename

    def get_basename(self) :
        return self.filename

    def get_extension(self) :
        return self.extension

class SffFile(DataFile) :
    def __init__(self, fname) :
        DataFile.__init__(self, fname, ".sff")

class State(object) :
    def __init__(self, states) :
        self.__counter = 0
        self.__num_states = states

    def inc(self) :
        self.__counter += 1
        self.__counter %= self.__num_states

    def get(self) :
        return self.__counter

    def __len__(self) :
        return self.__num_states

class ParseError(Exception) :
    pass

class FastqFile(DataFile) :
    SEQID = 0
    SEQ = 1
    QUALID = 2
    QUAL = 3

    def __init__(self, fname) :
        super(FastqFile, self).__init__(fname, ".fastq")
        self._filehandle = None #open(self.get_filename())
        self._state = FastqFile.SEQID
        self._linenum = 0

        self._validators = {
                FastqFile.SEQID  : self.__validate_seqid,
                FastqFile.SEQ    : self.__validate_sequence,
                FastqFile.QUALID : self.__validate_qualid,
                FastqFile.QUAL   : self.__validate_qualities
            }

        self._current = {  
                FastqFile.SEQID  : None,
                FastqFile.SEQ    : None,
                FastqFile.QUALID : None,
                FastqFile.QUAL   : None 
            }

    def __validate_seqid(self, s) :
        if not (s.startswith('@') or s.startswith('>')) :
            raise ParseError("%s : expected line %d to start with a @ or > (started with %s)" % \
                    (self.get_filename(), self._linenum, s[0]))

    def __validate_sequence(self, s) :
        uniq = set(s)

        for i in uniq :
            if i not in IUPAC.codes :
                raise ParseError("%s : line %d contained an invalid UIPAC code (%s)" % \
                        (self.get_filename(), self._linenum, i))

    def __validate_qualid(self, s) :
        if not s.startswith('+') :
            raise ParseError("%s : expected line %d to start with a +" % \
                    (self.get_filename(), self._linenum))

    def __validate_qualities(self, s) :
        uniq = set(s)

        for i in uniq :
            try :
                Sequence.quality_to_int(i)
        
            except ValueError, ve :
                raise ParseError("%s : line %d contained an invalid quality value (%s)" % \
                        (self.get_filename(), self._linenum, i))

    def __validate(self, s, state) :
        self._validators[state](s)

    def __iter__(self) :
        return self

    def next(self) :
        return self.read()

    def open(self) :
        if self._filehandle :
            self._filehandle.close()

        self._filehandle = open(self.get_filename())

    def close(self) :
        if self._filehandle :
            self._filehandle.close()

    def seq(self) :
        seqid = self._current[FastqFile.SEQID]
        duplicates = 1

        if "NumDuplicates" in self._current[FastqFile.SEQID] :
            mat = re.match(">(\S+)\ NumDuplicates=(\d+)$", self._current[FastqFile.SEQID])

            if not mat :
                raise DataFileError("'%s' is a malformed sequence id" % self._current[FastqFile.SEQID])

            seqid = mat.group(1)
            duplicates = int(mat.group(2))

        tmp = Sequence(self._current[FastqFile.SEQ], 
                None if self._current[FastqFile.QUAL] == "" else self._current[FastqFile.QUAL])
        
        # hack, maybe make more documented
        tmp.id = seqid[1:]
        tmp.duplicates = duplicates

        return tmp

    def read(self) :
        for line in self._filehandle :
            line = line.strip()

            self._linenum += 1

            if line == "" :
                continue

            # both '@' and '>' are legitimate quality scores
            # but '+' is a genuine delimiter
            # edit: unfortunately so is '+', so i also need to be in the SEQ state
            if line.startswith('+') and (self._state == FastqFile.SEQ) :
                self._current[FastqFile.QUALID] = line
                self._state = FastqFile.QUAL
                continue

            if self._state == FastqFile.SEQID :
                self.__validate(line, self._state)

                self._current[FastqFile.SEQID] = line
                self._current[FastqFile.SEQ] = ""
                self._current[FastqFile.QUALID] = ""
                self._current[FastqFile.QUAL] = ""

                self._state = FastqFile.SEQ

            elif self._state == FastqFile.SEQ :
                # if we are reading a fasta file
                if line.startswith('>') :
                    tmp = self.seq()
                    self._current[FastqFile.SEQID] = line
                    self._current[FastqFile.SEQ] = ""
                    return tmp

                self._current[FastqFile.SEQ] += line
                # note: sequence can be across an arbitrary number of lines, so we have
                # to stay in state SEQ 

            elif self._state == FastqFile.QUAL :
                self.__validate(line, self._state)

                self._current[FastqFile.QUAL] += line

                if len(self._current[FastqFile.SEQ]) == len(self._current[FastqFile.QUAL]) :
                    self._state = FastqFile.SEQID
                    return self.seq()

        if self._current[FastqFile.SEQID] is not None and self._current[FastqFile.SEQID].startswith('>') :
            tmp = self.seq()
            self._current[FastqFile.SEQID] = ""
            return tmp

        raise StopIteration

