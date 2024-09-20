#!/usr/bin/env python3

from struct import unpack, pack, calcsize, unpack_from
import numpy as np
from math import cos, sqrt
from zlib import crc32
from string import printable
from os.path import isfile
from glob import glob

class JOBARC(object):
    MAX_RECORDS = 1000
    RECORD_SIZE = 128

    class Record(object):
        def __init__(self, label, start, size):
            self.label = label.decode("utf-8")
            self.start = start
            self.size = size

    @staticmethod
    def parse_occ(s):
        occ = [[int(o) for o in half.split('-')] for half in s.split('/', 1)]
        if len(occ) == 1:
            occ.append(occ[0])
        return occ

    @staticmethod
    def parse_orbs(s):
        orbs = []
        for o in s.split('-'):
            if '>' in o:
                frm, to = o.split('>')
                orbs += range(int(frm), int(to)+1)
            else:
                orbs.append(int(o))
        return orbs

    def __init__(self):
        self.records = []

        fd = open('JAINDX', 'rb')
        data = fd.read()
        fd.close()

        if len(data) == 8000 + 2001*4 + 2*4:
            self.isize = 4
            self.msize = 4
        elif len(data) == 8000 + 2001*8 + 2*4:
            self.isize = 8
            self.msize = 4
        elif len(data) == 8000 + 2001*4 + 2*8:
            self.isize = 4
            self.msize = 8
        elif len(data) == 8000 + 2001*8 + 2*8:
            self.isize = 8
            self.msize = 8
        else:
            raise RuntimeError("JAINDX file is corrupt")

        self.itype = 'q' if self.isize == 8 else 'i'
        self.mtype = 'q' if self.msize == 8 else 'i'

        off = self.msize
        labels = [data[off+i*8:off+i*8+8] for i in range(1000)]
        off += 8000
        start = unpack('=1000'+self.itype, data[off:off+self.isize*1000])
        off += self.isize*1000
        size = unpack('=1000'+self.itype, data[off:off+self.isize*1000])

        for i in range(1000):
            if labels[i] != 'OPENSLOT':
                self.records.append(self.Record(labels[i], start[i], size[i]))

        self.fd = open('JOBARC', 'rb')

        self.iflags = np.append(self.getrec("IFLAGS", self.itype, 100),
                                self.getrec("IFLAGS2", self.itype, 500))

    def __del__(self):
        self.fd.close()

    def record(self, label):
        for rec in self.records:
            if rec.label == label:
                return rec

    def getrec(self, label, type, num=None):
        label = label[:8].ljust(8)
        rec = self.record(label)

        if not rec:
            raise RuntimeError("Record " + label + " not found")

        real_num = num if num else 1
        format = '=' + str(real_num) + type
        bytes = calcsize(format)

        if bytes > rec.size*self.isize:
            raise RuntimeError("Attempting to read past end of record %s (%d vs. %d)" % (label, bytes, rec.size*self.isize))

        self.fd.seek((rec.start-1)*self.isize, 0)

        if not num:
            return unpack(format, self.fd.read(bytes))[0]
        else:
            type_map = {'d': 'float64', 'f': 'float32', 'i': 'int32', 'q': 'int64'}
            return np.fromfile(self.fd, type_map[type], num)

    def putrec(self, label, data):
        label = label[:8].ljust(8)
        rec = self.record(label)

        def typestr(val):
            if isinstance(val, int):
                return self.itype
            elif isinstance(val, float):
                return 'd'
            else:
                raise TypeError("Wrong type for list data")

        if isinstance(data, (str, unicode)):
            data = data.encode()
        elif isinstance(data, (list, tuple)):
            data = pack('=' + str(len(data)) + typestr(data[0]), *data)
        else:
            data = pack('=' + typestr(data), data)

        if not rec:
            if len(self.records) == self.MAX_RECORDS:
                raise RuntimeError("Maximum number of records reached")

            start = max([rec.start + rec.size for rec in self.records])
            size = (len(data) + self.isize - 1)//self.isize
            rec = self.Record(label, start, size)
            self.records.append(rec)

        if len(data) > rec.size*self.isize:
            raise RuntimeError("Attempting to write past end of record %s (%d vs. %d)" % (label, len(data), rec.size*self.isize))

        self.fd.seek((rec.start-1)*self.isize, 0)
        self.fd.write(data)

    def recsize(self, label):
        label = label[:8].ljust(8)
        rec = self.record(label)
        return rec.size*self.isize if rec else -1

    def recexists(self, label):
        return self.recsize(label) != -1
        
    def getflag(self, flag):
        return self.iflags[flag-1]
        
class MOINTS(object):
    filenames = ("MOINTS", "GAMLAM", "MOABCD", "DERINT", "DERGAM")
        
    def __init__(self, jobarc):
        self.moio = jobarc.getrec("MOIOVEC", jobarc.itype, 5000)
        self.moiowd = jobarc.getrec("MOIOWRD", jobarc.itype, 5000)
        self.moiofl = jobarc.getrec("MOIOFIL", jobarc.itype, 5000)
        self.moiods = jobarc.getrec("MOIODIS", jobarc.itype, 5000)
        self.moiosz = jobarc.getrec("MOIOSIZ", jobarc.itype, 5000)
        self.isytyp = jobarc.getrec("ISYMTYP", jobarc.itype, 1000)
        self.record_size = jobarc.getflag(37)
        self.isize = jobarc.isize
        
    def read(self, list, off=0):
        j = (list-1)*10+off
        
        sizel = self.moiosz[j]
        sizer = self.moiods[j]
        n = sizel*sizer

        off = ((self.moio[j]-1)*self.record_size + self.moiowd[j]-1)*self.isize
        
        with open(self.filenames[(list-1)//100], 'rb') as fd:
            fd.seek(off, 0)
            return np.fromfile(fd, 'float64', n)

if __name__ == '__main__':
    jobarc = JOBARC()
    moints = MOINTS(jobarc)

    nirrep = jobarc.getrec("COMPNIRR", jobarc.itype);
    assert(nirrep == 1)
    nocc = jobarc.getrec("SYMPOPOA", jobarc.itype)
    nvrt = jobarc.getrec("SYMPOPVA", jobarc.itype)
    
    fi = np.zeros(nocc*nocc)
    fa = np.zeros(nvrt*nvrt)
    
    evals = jobarc.getrec("SCFEVALA", 'd', nocc+nvrt)
    energy = jobarc.getrec("TOTENERG", 'd') - \
             jobarc.getrec("SCFENEG ", 'd')
    
    for i in range(nocc):
        fi[i+i*nocc] = evals[i]
    
    for i in range(nvrt):
        fa[i+i*nvrt] = evals[nocc+i]

    with open('mp2_energy', 'w') as fd:
        fd.write('%.15f\n' % (energy))

    with open('fi.dat', 'wb') as fd:
        fd.write(pack('=i', nocc))
        fi.tofile(fd)

    with open('fa.dat', 'wb') as fd:
        fd.write(pack('=i', nvrt))
        fa.tofile(fd)
        
    v = moints.read(16)
    assert(len(v) == nvrt*nvrt*nocc*nocc)

    with open('v.dat', 'wb') as fd:
        v.tofile(fd)

    t = moints.read(46)
    assert(len(t) == nvrt*nvrt*nocc*nocc)

    with open('t.dat', 'wb') as fd:
        t.tofile(fd)

