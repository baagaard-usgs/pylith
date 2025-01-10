#!/usr/bin/env nemesis

# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

import numpy
import h5py
import matplotlib.pyplot as plt
from axialtraction_maxwell_soln import AnalyticalSoln

from pythia.pyre.units.time import year

# Input files.
lineDefs = ['r+', 'g+', 'b+', 'c+', 'm+', 'y+']
lineDefAnl = 'k-'
legendAnl = 'Analytical'
basePrefix = 'axialtraction_powerlaw_n1_dt'

    
def getVars(fileName):
    """
    Get time, sxx, eyy, and dispy from given file.
    """
    h5 = h5py.File(fileName, 'r')
    time = h5['time'][:].flatten()
    timeYears = time / year.value
    locs = h5['geometry/vertices'][:]
    stress = h5['vertex_fields/cauchy_stress'][:]
    strain = h5['vertex_fields/cauchy_strain'][:]
    disp = h5['vertex_fields/displacement'][:]
    sxx = stress[:,0,0]
    eyy = strain[:,0,1]
    dispy = disp[:,0,1]

    h5.close()

    return (timeYears, sxx, eyy, dispy, locs)


def scanLogfile(fileName):
    """
    Get Jacobian info from log file.
    """
    matchStr = '||J - Jfd||_F/||J||_F ='
    endStr = ','
    lastLine = None
    val = None
    try:
        f = open(fileName, 'r')
        lines = reversed(f.readlines())
        for line in lines:
            if (matchStr in line):
                lastLine = line
                break
        if (lastLine):
            val = float((lastLine.split(matchStr))[1].split(endStr)[0])
    except FileNotFoundError:
        pass

    return val

    
def run(axs, colNum, stepSizes, useJacobian):
    """
    Create subplots and loop over simulations.
    """
    numSims = len(stepSizes)
    jacobianDiff = numpy.zeros(numSims, dtype=numpy.float64)
    jacobianInfo = False

    for simNum in range(numSims):
        dt = stepSizes[simNum]
        dtStr = f"{dt:04.2f}"
        baseName = basePrefix + dtStr
        h5File = 'output/' + baseName + '-viscomat.h5'
        logFile = baseName + '.log'
        (timeYears, sxx, eyy, dispy, locs) = getVars(h5File)
        timeSecs = timeYears*year.value
        jacobianDiff[simNum] = scanLogfile(logFile)
        if (jacobianDiff[simNum]):
            jacobianInfo = True
        if (simNum == 0):
            (dispyAnl, stressAnl, devStressAnl, strainAnl, devStrainAnl,
             maxwellVisStrain, powerLawVisStrain) = AnalyticalSoln(timeSecs, locs[0,:].reshape(1,2))
            axs[0][colNum].plot(timeYears, stressAnl[:,0], lineDefAnl, label=legendAnl)
            axs[1][colNum].plot(timeYears, strainAnl[:,1], lineDefAnl, label=legendAnl)
            axs[2][colNum].plot(timeYears, dispyAnl, lineDefAnl, label=legendAnl)
        axs[0][colNum].plot(timeYears, sxx, lineDefs[simNum], label="dt={0}".format(dtStr))
        axs[1][colNum].plot(timeYears, eyy, lineDefs[simNum], label="dt={0}".format(dtStr))
        axs[2][colNum].plot(timeYears, dispy, lineDefs[simNum], label="dt={0}".format(dtStr))
        
    axs[0][colNum].set_xlabel('Time (years)')
    axs[0][colNum].set_ylabel('Stress_xx (Pa)')
    axs[0][colNum].legend(loc="upper right")
    axs[0][colNum].set_title('Axial traction n=1')
    axs[1][colNum].set_xlabel('Time (years)')
    axs[1][colNum].set_ylabel('Strain_yy')
    axs[1][colNum].legend(loc="upper right")
    axs[2][colNum].set_xlabel('Time (years)')
    axs[2][colNum].set_ylabel('Displacement_y (m)')
    axs[2][colNum].legend(loc="upper right")
    jacobianPlot = useJacobian and jacobianInfo
    if (jacobianPlot):
        axs[3][colNum].loglog(stepSizes, jacobianDiff, 'k+-')
        axs[3][colNum].set_xlabel('Time step size (years)')
        axs[3][colNum].set_ylabel('Jacobian difference')


# End of file
