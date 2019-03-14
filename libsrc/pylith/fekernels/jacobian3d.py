#!/usr/bin/env python
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2015 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#
# Initial attempt to compute plane strain Jacobian matrices symbolically.
# PREREQUISITES:  sympy
# ----------------------------------------------------------------------
#
# import pdb
# pdb.set_trace()
import sympy
import sympy.tensor
import sympy.tensor.array
from itertools import product
# ----------------------------------------------------------------------
ndim = 3
numComps = 3
ndimRange = range(ndim)
numCompsRange = range(numComps)

# Constants.
zero = sympy.sympify(0)
one = sympy.sympify(1)
two = sympy.sympify(2)
three = sympy.sympify(3)

# Define basis and displacement vector.
from sympy.abc import x, y, z
u1, u2, u3 = sympy.symbols('u1 u2 u3', type="Function")
X = [x, y, z]
U = [u1(x,y,z), u2(x,y,z), u3(x,y,z)]

# Deformation gradient, transpose, and strain tensor.
defGrad = sympy.derive_by_array(U, X)
defGradTranspose = defGrad.transpose()
strain = ((defGrad + defGradTranspose)/two).tomatrix()

# Define volumetric strain and deviatoric strain.
volStrain = sympy.tensorcontraction(strain, (0, 1))
volStrainArr = volStrain * sympy.eye(ndim)
devStrain = strain - volStrainArr/three

# Define displacements and strains for previous time step.
un1, un2, un3 = sympy.symbols('un1 un2 un3', type="Function")
Un = [un1(x,y,z), un2(x,y,z), un3(x,y,z)]
defGradN = sympy.derive_by_array(Un, X)
defGradNTranspose = defGradN.transpose()
strainN = ((defGradN + defGradNTranspose)/two).tomatrix()
meanStrainN = sympy.tensorcontraction(strainN, (0, 1))/three
meanStrainNArr = meanStrainN * sympy.eye(ndim)
devStrainN = strainN - meanStrainNArr

# Put in dummy tensor for things that don't depend on U.
aa, ab, ac, ba, bb, bc, ca, cb, cc = sympy.symbols('aa ab ac ba bb bc ca cb cc')
dummyTensor = sympy.Matrix([[aa, ab, ac],
                            [ba, bb, bc],
                            [ca, cb, cc]])

# ----------------------------------------------------------------------
def divideTensors(t1, t2):
  """
  Function to divide each element of t1 by the corresponding element of t2.
  """
  from sympy import MutableDenseNDimArray
  tRank = t1.rank()
  tSize = t1.shape[0]
  tLen = tSize**tRank
  tReturn = MutableDenseNDimArray(range(tLen), t2.shape)
  for elem in range(tLen):
    tReturn[elem] = t1[elem]/t2[elem]

  return tReturn


def innerProd(mat1, mat2):
  """
  Function to compute the scalar inner product of two matrices.
  I am sure there is a much easier way to do this.
  """
  test1 = isinstance(mat1, sympy.MutableDenseMatrix)
  test2 = isinstance(mat2, sympy.MutableDenseMatrix)
  m1 = mat1.copy()
  m2 = mat2.copy()
  if (not test1):
    m1 = mat1.tomatrix()
  if (not test2):
    m2 = mat2.tomatrix()
  matMult = sympy.matrix_multiply_elementwise(m1, m2)
  rowSum = matMult * sympy.ones(matMult.shape[1], 1)
  colSum = sympy.ones(1, rowSum.shape[0]) * rowSum
  scalarProd = colSum[0]

  return scalarProd


def writeJacobianUniqueVals(f, jacobian):
  """
  Function to write unique values and assign them to variables.
  """

  # Unique entries in Jacobian, excluding 0.
  uniqueVals = list(set(jacobian))
  uniqueVals.remove(0)
  numUniqueVals = len(uniqueVals)
  uniqueValNames = numUniqueVals * [None]
  usedVals = numUniqueVals * [None]
  
  f.write("/* Unique components of Jacobian. */\n")
  outFmt = "const PylithReal %s = %s;\n"

  # Loop over Jacobian components in original PyLith order.
  ui = 0
  for i, j, k, l in product(numCompsRange, ndimRange, numCompsRange, ndimRange):
    ii = i + 1
    jj = j + 1
    kk = k + 1
    ll = l + 1
    if (jacobian[i,j,k,l] in uniqueVals and jacobian[i,j,k,l] not in usedVals):
      testInd = uniqueVals.index(jacobian[i,j,k,l])
      comp = "C" + repr(ii) + repr(jj) + repr(kk) + repr(ll)
      f.write(outFmt % (comp, jacobian[i,j,k,l]))
      uniqueValNames[testInd] = comp
      usedVals[ui] = jacobian[i,j,k,l]
      ui += 1
    if (ui == numUniqueVals):
      break

  return (uniqueVals, uniqueValNames)


def writeJacobianComments(f, jacobian):
  """
  Function to write correspondence between PETSc and PyLith Jacobian values.
  """

  f.write("/* j(f,g,df,dg) = C(f,df,g,dg)\n\n")
  outFmt = "%d:  %s = %s = %s\n"
  
  # Loop over Jacobian components in new order.
  ui = 0
  for i, k, j, l in product(numCompsRange, numCompsRange, ndimRange, ndimRange):
    ii = i + 1
    jj = j + 1
    kk = k + 1
    ll = l + 1
    pyComp = "C" + repr(ii) + repr(jj) + repr(kk) + repr(ll)
    peComp = "j" + repr(i) + repr(k) + repr(j) + repr(l)
    f.write(outFmt % (ui, peComp, pyComp, jacobian[i,j,k,l]))
    ui += 1
            
  f.write("*/\n\n")

  return


def writeJacobianNonzero(f, jacobian, uniqueVals, uniqueValNames):
  """
  Function to write nonzero Jacobian entries using predefined value names.
  """

  f.write("/* Nonzero Jacobian entries. */\n")
  
  outFmt = "Jg3[%d] -=  %s; /* %s */\n"
  
  # Loop over Jacobian components in new order.
  ui = 0
  for i, k, j, l in product(numCompsRange, numCompsRange, ndimRange, ndimRange):
    ii = i + 1
    jj = j + 1
    kk = k + 1
    ll = l + 1
    peComp = "j" + repr(i) + repr(k) + repr(j) + repr(l)
    if (jacobian[i,j,k,l] != 0):
      ind = uniqueVals.index(jacobian[i,j,k,l])
      f.write(outFmt % (ui, uniqueValNames[ind], peComp))

    ui += 1

  return


def writeJacobianInfo(fileName, jacobian):
  """
  Function to write info about Jacobian.
  """
  f = open(fileName, 'w')

  (uniqueVals, uniqueValNames) = writeJacobianUniqueVals(f, jacobian)
  writeJacobianComments(f, jacobian)
  writeJacobianNonzero(f, jacobian, uniqueVals, uniqueValNames)
  f.close()

  return
                  

# ----------------------------------------------------------------------
# Elastic isotropic stress.
fileName = 'elasticity-elas_iso3d.txt'
(lambdaModulus, shearModulus,
 bulkModulus) = sympy.symbols('lambdaModulus shearModulus bulkModulus')
stress = lambdaModulus * volStrainArr + two * shearModulus * strain
meanStress = sympy.tensorcontraction(stress, (0, 1))/three
meanStressArr = meanStress * sympy.eye(ndim)
devStress = stress - meanStressArr
jacobian = sympy.derive_by_array(stress, defGrad)
writeJacobianInfo(fileName, jacobian)

# ----------------------------------------------------------------------
# Maxwell viscoelastic.
fileName = 'elasticity-max_iso3d.txt'
deltaT, tauM = sympy.symbols('deltaT tauM')
expFac = sympy.exp(-deltaT/tauM)
dq = sympy.symbols('dq')
delHArr = dq * (devStrain - devStrainN)
# Dummy tensor represents viscous strain from previous time step.
hMArr = expFac * dummyTensor + delHArr
meanStress = bulkModulus * volStrainArr
devStress = two * shearModulus * hMArr
stress = meanStress + devStress
jacobian = sympy.derive_by_array(stress, defGrad)
writeJacobianInfo(fileName, jacobian)

# ----------------------------------------------------------------------
# Generalized Maxwell viscoelastic.
fileName = 'elasticity-genmax_iso3d.txt'
(tauM1, tauM2, tauM3, shearModulusRatio_1, shearModulusRatio_2,
 shearModulusRatio_3) = sympy.symbols(
  'tauM1 tauM2 tauM3 shearModulusRatio_1 shearModulusRatio_2 shearModulusRatio_3')
shearModulusRatio_0 = sympy.symbols('shearModulusRatio_0')
expFac1 = sympy.exp(-deltaT/tauM1)
expFac2 = sympy.exp(-deltaT/tauM2)
expFac3 = sympy.exp(-deltaT/tauM3)
dq_1, dq_2, dq_3 = sympy.symbols('dq_1 dq_2 dq_3')
delHArr1 = dq_1 * (devStrain - devStrainN)
delHArr2 = dq_2 * (devStrain - devStrainN)
delHArr3 = dq_3 * (devStrain - devStrainN)
# Dummy tensors represent viscous strain from previous time step.
hMArr1 = expFac1 * dummyTensor + delHArr1
hMArr2 = expFac2 * dummyTensor + delHArr2
hMArr3 = expFac3 * dummyTensor + delHArr3
meanStress = bulkModulus * volStrainArr
devStress = two * shearModulus * (shearModulusRatio_0 * devStrain + \
                                  shearModulusRatio_1 * hMArr1 + \
                                  shearModulusRatio_2 * hMArr2 + \
                                  shearModulusRatio_3 * hMArr3)
stress = meanStress + devStress
jacobian = sympy.derive_by_array(stress, defGrad)
writeJacobianInfo(fileName, jacobian)

# ----------------------------------------------------------------------
# Power-law viscoelastic.
fileName = 'elasticity-powerlaw_iso3d.txt'
aT, n, alpha = sympy.symbols('aT n alpha')
s11, s12, s13, s22, s23, s33 = sympy.symbols('s11 s12 s13 s22 s23 s33')
s11T, s12T, s13T, s22T, s23T, s33T = sympy.symbols(
  's11T s12T s13T s22T s23T s33T')
s11FTau, s12FTau, s13FTau, s22FTau, s23FTau, s33FTau = sympy.symbols(
  's11FTau s12FTau s13FTau s22FTau s23FTau s33FTau')
j2FTplusDt, j2Ft, j2FTau, gammaFTau = sympy.symbols(
  'j2FTplusDt j2Ft j2FTau gammaFTau')
meanStress = bulkModulus * volStrainArr
devStress = sympy.Matrix([[s11, s12, s13],
                          [s12, s22, s23],
                          [s13, s23, s33]])
devStressT = sympy.Matrix([[s11T, s12T, s13T],
                           [s12T, s22T, s23T],
                           [s13T, s23T, s33T]])
devStressTau = alpha*devStress + (one - alpha)*devStressT
j2TplusDt = sympy.sqrt(innerProd(devStress, devStress))
j2T = sympy.sqrt(innerProd(devStressT, devStressT))
j2Tau = sympy.sqrt(innerProd(devStressTau, devStressTau))
aE = one/(two*shearModulus)
gammaTau = aT*(j2Tau)**(n-one)
F = aE*devStress + devStressTau*gammaTau*deltaT - devStrain

dFdStress = sympy.derive_by_array(F, devStress)
dFdStressSimp = dFdStress.subs([(gammaTau, gammaFTau),
                                (j2Tau, j2FTau),
                                (devStressTau[0], s11FTau),
                                (devStressTau[1], s12FTau),
                                (devStressTau[2], s13FTau),
                                (devStressTau[4], s22FTau),
                                (devStressTau[5], s23FTau),
                                (devStressTau[8], s33FTau)])

dFdStrain = -sympy.derive_by_array(F, defGrad)

jacobianDev = divideTensors(dFdStrain, dFdStressSimp)
jacobianVol = sympy.derive_by_array(meanStress, defGrad)
jacobian = jacobianDev + jacobianVol
writeJacobianInfo(fileName, jacobian)
