import os
import re

import numpy as np
import cclib


class XtbLog:
    def __init__(self, file):
        # default values for thermochemical calculations
        if '.log' not in file:
            raise TypeError('A xtb .log file must be provided')

        self.file = file
        self.name = os.path.basename(file)

        self.GetTermination()
        if not self.termination:
            pass
            #self.GetError()
        else:
            self.GetFreq()
            self.GetE()

    def GetTermination(self):
        with open(self.file) as fh:
            for line in (fh):
                if line.find("normal termination") > -1:
                    self.termination = True
                    return True
            self.termination = False

    def GetFreq(self):
        with open(self.file) as fh:
            txt = fh.readlines()

        txt = [x.strip() for x in txt]
        for i, line in enumerate(txt):
            if line.find('Frequency Printout') > -1:
                txt = txt[i + 3:]
                break

        waveNums = []
        for i, line in enumerate(txt):
            if line.find('reduced masses') > -1:
                txt = txt[i + 1:]
                break
            m = re.findall('\s+(-?\d+\.\d+)', line)
            if m:
                for match in m:
                    waveNums.append(float(match.strip()))

        for i, line in enumerate(txt):
            if line.find('IR intensities') > -1:
                txt = txt[i + 1:]
                break

        intensities = []
        for i, line in enumerate(txt):
            if line.find('Raman intensities') > -1:
                txt = txt[i + 1:]
                break
            m = re.findall('\d+:\s+(\d+\.\d+)', line)
            if m:
                for match in m:
                    intensities.append(float(match))

        waveNums, intensities = list(zip(*[(w, i) for w, i in zip(waveNums, intensities) if w != 0]))

        if waveNums and intensities and len(waveNums) == len(intensities):
            self.wavenum = waveNums
            self.ir_intensities = intensities

    def GetE(self):
        with open(self.file) as fh:
            txt = fh.readlines()

        txt = [x.strip() for x in txt]
        for i, line in enumerate(txt):
            m = re.search('TOTAL ENERGY\s+(-?\d+\.\d+)', line)
            if m:
                self.E = m[1]
                continue
            m = re.search('TOTAL ENTHALPY\s+(-?\d+\.\d+)', line)
            if m:
                self.H = m[1]
                continue
            m = re.search('TOTAL FREE ENERGY\s+(-?\d+\.\d+)', line)
            if m:
                self.G = float(m[1])


class G16Log:
    def __init__(self, file):
        # default values for thermochemical calculations
        if '.log' not in file:
            raise TypeError('A g16 .log file must be provided')

        self.file = file
        self.name = os.path.basename(file)

        self.GetTermination()
        if not self.termination:
            self.GetError()
        else:
            self.cclib()
            self.GetNMR()

    def GetTermination(self):
        with open(self.file) as fh:
            txt = fh.readlines()
            txt = [x.strip() for x in txt]
            if txt[-1].find("Normal termination") > -1:
                self.termination = True
                return True
            self.termination = False

    def GetError(self):
        with open(self.file) as fh:
            for line in (fh):
                if line.find("Error termination") > -1:
                    self.error = line
                    return True
            self.error = None

    def cclib(self):
        data = cclib.io.ccread(self.file)
        self.coords = data.atomcoords[-1]
        self.natom = data.natom
        self.scf = data.scfenergies[-1]

    def GetNMR(self):
        NMR = []
        with open(self.file, 'r') as fh:
            for line in fh:
                if len(NMR) == self.natom: break
                m = re.search('Isotropic\s*=\s*(-?\d+\.\d+)', line)
                if not m: continue
                NMR.append(float(m.group(1)))
        self.NMR = np.array(NMR)