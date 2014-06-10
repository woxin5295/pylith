cell = "hex8"
testCase = "open"

import numpy

from numpy import *
from numpy.linalg import inv

numpy.set_printoptions(precision=12)

lengthScale = 1.0e+3

# ----------------------------------------------------------------------
def printdata(data):
    """
    Print data as C array.
    """
    (nrows, ncols) = data.shape
    style = " %18.10e,"*ncols
    for row in xrange(nrows):
        print (style % tuple(data[row,:]))
    return


# ----------------------------------------------------------------------
def globalToFault(v, R):
    """
    Convert vector from global coordinate system to fault coordinate system.
    """
    (m,ndof) = v.shape

    vF = numpy.dot(R, v.reshape(m*ndof,1))
    return vF.reshape((m, ndof))


# ----------------------------------------------------------------------
def faultToGlobal(v, R):
    """
    Convert vector from fault coordinate system to global coordinate system.
    """
    (m,ndof) = v.shape

    vG = numpy.dot(R.transpose(), v.reshape(m*ndof,1))
    return vG.reshape((m, ndof))


# ----------------------------------------------------------------------
if cell in ["tri3", "tri3d", "quad4"]:
    if cell == "tri3":
        indexL = numpy.arange(6,8)
        indexN = numpy.arange(1,3)
        indexP = numpy.arange(4,6)
        n = 8
        m = 2
        DOF = 2

        A = numpy.array([[1.0, 0.0,],
                         [0.0, 1.0,],]);
        C = numpy.array([[0.0, +1.0, 0.0, 0.0,],
                         [+1.0, 0.0, 0.0, 0.0,],
                         [0.0, 0.0, 0.0, +1.0,],
                         [0.0, 0.0, +1.0, 0.0,],]);
    
        fieldT = numpy.array([[ 8.1, 9.1,],
                              [ 8.2, 9.2,],
                              [ 8.3, 9.3,],
                              [ 8.4, 9.4,],
                              [ 8.2, 9.2,],
                              [ 8.3, 9.3,],
                              [-8.6, 9.6,],
                              [-8.8, 9.8,],])

        if testCase == "stick":
            fieldTIncr = numpy.array([[ 1.1, 2.1,],
                                      [ 1.2, 2.2,],
                                      [ 1.3, 2.3,],
                                      [ 1.4, 2.4,],
                                      [ 1.2, 2.2,],
                                      [ 1.3, 2.3,],
                                      [-21.6, 2.6,],
                                      [-21.8, 2.8,],])            
        elif testCase == "slip":
            fieldTIncr = numpy.array([[ 1.1, 2.1,],
                                      [ 1.2, 2.2,],
                                      [ 1.3, 2.3,],
                                      [ 1.4, 2.4,],
                                      [ 1.2, 2.2,],
                                      [ 1.3, 2.3,],
                                      [-1.6, 2.6,],
                                      [-1.8, 2.5,],])            
        elif testCase == "open":
            fieldTIncr = numpy.array([[ 9.1, 7.1,],
                                      [ 9.2, 7.2,],
                                      [ 9.3, 7.3,],
                                      [ 9.4, 7.4,],
                                      [ 9.5, 7.5,],
                                      [ 9.6, 7.6,],
                                      [ +10.6, -10.6,],
                                      [ +10.8, -10.8,],])            


    elif cell == "tri3d":
        indexL = numpy.array([9, 10, 11])
        indexN = numpy.array([1, 2, 4])
        indexP = numpy.array([6, 7, 8])
        n = 12
        m = 3
        DOF = 2

        A = numpy.array([[2.0, 0.0, 0.0],
                         [0.0, 1.0, 0.0],
                         [0.0, 0.0, 1.0],])
        C = numpy.array([[-0.70710678118654757, +0.70710678118654757, 0.0, 0.0, 0.0, 0.0,],
                         [+0.70710678118654757, +0.70710678118654757, 0.0, 0.0, 0.0, 0.0,],
                         [0.0, 0.0, 0.0, +1.0, 0.0, 0.0,],
                         [0.0, 0.0, +1.0, 0.0, 0.0, 0.0,],
                         [0.0, 0.0, 0.0, 0.0, -1.0, 0.0,],
                         [0.0, 0.0, 0.0, 0.0, 0.0, +1.0,],])
    
        fieldT = numpy.array([[ 6.1, 8.1,],
                              [ 6.2, 8.2,],
                              [ 6.3, 8.3,],
                              [ 6.4, 8.4,],
                              [ 6.5, 8.5,],
                              [ 6.6, 8.6,],
                              [ 6.2, 8.2,],
                              [ 6.3, 8.3,],
                              [ 6.5, 8.5,],
                              [-3.8,-4.8,],
                              [-3.0, 4.0,],
                              [ 3.2,-4.2,],])

        if testCase == "stick":
            fieldTIncr = numpy.array([[ 1.1, 2.1,],
                                      [ 1.2, 2.2,],
                                      [ 1.3, 2.3,],
                                      [ 1.4, 2.4,],
                                      [ 1.5, 2.5,],
                                      [ 1.6, 2.6,],
                                      [ 1.2, 2.2,],
                                      [ 1.3, 2.3,],
                                      [ 1.5, 2.5,],
                                      [-21.8,-22.8,],
                                      [-21.0, 2.0,],
                                      [ 2.2,-22.2,],])            
        elif testCase == "slip":
            fieldTIncr = numpy.array([[ 1.1, 2.1,],
                                      [ 1.2, 2.2,],
                                      [ 1.3, 2.3,],
                                      [ 1.4, 2.4,],
                                      [ 1.5, 2.5,],
                                      [ 1.6, 2.6,],
                                      [ 1.2, 2.2,],
                                      [ 1.3, 2.3,],
                                      [ 1.5, 2.5,],
                                      [-1.8,+3.6,],
                                      [-1.0, 1.1,],
                                      [ 1.7,-1.2,],])            
        elif testCase == "open":
            fieldTIncr = numpy.array([[ 1.1, 2.1,],
                                      [ 1.2, 2.2,],
                                      [ 1.3, 2.3,],
                                      [ 1.4, 2.4,],
                                      [ 1.5, 2.5,],
                                      [ 1.6, 2.6,],
                                      [ 1.2, 2.2,],
                                      [ 1.3, 2.3,],
                                      [ 1.5, 2.5,],
                                      [+11.8, 11.8,],
                                      [+10.0, 0.1,],
                                      [ 1.2, +10.2,],])            


    elif cell == "quad4":
        indexL = numpy.arange(8,10)
        indexN = numpy.arange(2,4)
        indexP = numpy.arange(6,8)
        n = 10
        m = 2
        DOF = 2

        A = numpy.array([[1.0, 0.0,],
                         [0.0, 1.0,],]);
        C = numpy.array([[0.0, +1.0, 0.0, 0.0,],
                         [+1.0, 0.0, 0.0, 0.0,],
                         [0.0, 0.0, 0.0, +1.0,],
                         [0.0, 0.0, +1.0, 0.0,],]);
    

        fieldT = numpy.array([[ 8.1, 9.1,],
                              [ 8.3, 9.3,],
                              [ 8.2, 9.2,],
                              [ 8.3, 9.3,],
                              [ 8.5, 9.5,],
                              [ 8.6, 9.6,],
                              [ 8.2, 9.6,],
                              [ 8.3, 9.8,],
                              [-8.6, 9.6,],
                              [-8.8, 9.8,],])

        if testCase == "stick":
            fieldTIncr = numpy.array([[ 1.1, 2.1,],
                                      [ 1.4, 2.4,],
                                      [ 1.2, 2.2,],
                                      [ 1.3, 2.3,],
                                      [ 1.5, 2.5,],
                                      [ 1.6, 2.6,],
                                      [ 1.2, 2.2,],
                                      [ 1.3, 2.3,],
                                      [-21.6, 2.6,],
                                      [-21.8, 2.5,],])
        elif testCase == "slip":
            fieldTIncr = numpy.array([[ 1.1, 2.1,],
                                      [ 1.2, 2.2,],
                                      [ 1.2, 2.2,],
                                      [ 1.3, 2.3,],
                                      [ 1.5, 2.5,],
                                      [ 1.6, 2.6,],
                                      [ 1.2, 2.2,],
                                      [ 1.3, 2.3,],
                                      [-1.6, 2.6,],
                                      [-1.8, 2.5,],])
        elif testCase == "open":
            fieldTIncr = numpy.array([[ 1.1, 2.1,],
                                      [ 1.2, 2.2,],
                                      [ 1.2, 2.2,],
                                      [ 1.3, 2.3,],
                                      [ 1.5, 2.5,],
                                      [ 1.6, 2.6,],
                                      [ 1.3, 2.7,],
                                      [ 1.4, 2.6,],
                                      [+10.6, -2.6,],
                                      [+10.8, -2.8,],])


    # ------------------------------------------------------------------
    A /= lengthScale

    fieldTpdt = (fieldT + fieldTIncr)
    fieldTpdtFault = globalToFault(fieldTpdt[indexL,:], C)

    tractionShearMag = numpy.abs(fieldTpdtFault[:,0])
    tractionNormal = fieldTpdtFault[:,1]

    print "fieldTpdt",fieldTpdt

    print "tractionShear",tractionShear
    print "tractionNormal",tractionNormal

    tractionRheologyFault = numpy.zeros((m, DOF))
    if testCase != "open":
        tractionRheologyShear = -0.6 * tractionNormal
        tractionRheologyNormal = tractionNormal

        tractionRheologyFault[:,0] = tractionRheologyShear * fieldTpdtFault[:,0]/tractionShearMag
        tractionRheologyFault[:,1] = tractionRheologyNormal
    print "tractionRheologyFault",tractionRheologyFault

    tractionRheology = faultToGlobal(tractionRheologyFault, C)
    tractionInternal = fieldTpdt[indexL,:]
    tractionResidual = tractionRheology - tractionInternal
    print "tractionRheology",tractionRheology
    print "tractionInternal",tractionInternal
    print "tractionResidual",tractionResidual

    residual = numpy.zeros(fieldT.shape)
    residual[indexN,:] = +numpy.dot(A, tractionResidual)
    residual[indexP,:] = -numpy.dot(A, tractionResidual)
    residual[indexL,:] = numpy.dot(A, tractionInternal) * (fieldTpdt[indexP,:] - fieldTpdt[indexN,:])

    print "residual \n",printdata(residual)


# ----------------------------------------------------------------------
elif cell in ["tet4", "hex8"]:

    if cell == "tet4":

        indexL = numpy.arange(8,11)
        indexN = numpy.arange(1,4)
        indexP = numpy.arange(5,8)
        n = 11
        m = 3
        DOF = 3

        A = numpy.array([[1.0/3.0,       0,        0,],
                         [      0, 1.0/3.0,        0,],
                         [      0,       0,  1.0/3.0,]])

        Cv = numpy.array([[ 0, +1, 0,],
                          [ 0, 0, +1,],
                          [ +1, 0, 0,],])
        Zv = numpy.zeros([3,3])
        C = numpy.vstack( (numpy.hstack((Cv, Zv, Zv)),
                           numpy.hstack((Zv, Cv, Zv)),
                           numpy.hstack((Zv, Zv, Cv)) ) )

        fieldT = numpy.array([[ 7.1, 8.1, 9.1,],
                              [ 7.2, 8.2, 9.2,],
                              [ 7.3, 8.3, 9.3,],
                              [ 7.4, 8.4, 9.4,],
                              [ 7.5, 8.5, 9.5,],
                              [ 7.2, 8.2, 9.2,],
                              [ 7.3, 8.3, 9.3,],
                              [ 7.4, 8.4, 9.4,],
                              [-7.7, 18.7, 19.7,],
                              [-7.9, 18.9, 19.9,],
                              [-7.1, 18.1, 19.1,],])

        if testCase == "stick":
            fieldTIncr = numpy.array([[ 1.1, 2.1, 3.1,],
                                      [ 1.2, 2.2, 3.2,],
                                      [ 1.3, 2.3, 3.3,],
                                      [ 1.4, 2.4, 3.4,],
                                      [ 1.5, 2.5, 3.5,],
                                      [ 1.2, 2.2, 3.2,],
                                      [ 1.3, 2.3, 3.3,],
                                      [ 1.4, 2.4, 3.4,],
                                      [-81.7, 2.7, 3.7,],
                                      [-81.9, 2.9, 3.9,],
                                      [-81.1, 2.1, 3.1,],])            
        elif testCase == "slip":
            fieldTIncr = numpy.array([[ 1.1, 2.1, 3.1,],
                                      [ 1.2, 2.2, 3.2,],
                                      [ 1.3, 2.3, 3.3,],
                                      [ 1.4, 2.4, 3.4,],
                                      [ 1.5, 2.5, 3.5,],
                                      [ 1.2, 2.5, 3.4,],
                                      [ 1.3, 2.4, 3.5,],
                                      [ 1.4, 2.6, 3.6,],
                                      [-4.7, 5.7, 6.7,],
                                      [-4.9, 5.9, 6.9,],
                                      [-4.1, 5.1, 6.1,],])            
        elif testCase == "open":
            fieldTIncr = numpy.array([[ 1.1, 2.1, 3.1,],
                                      [ 1.2, 2.2, 3.2,],
                                      [ 1.3, 2.3, 3.3,],
                                      [ 1.4, 2.4, 3.4,],
                                      [ 1.5, 2.5, 3.5,],
                                      [ 1.5, 3.2, 4.2,],
                                      [ 1.4, 3.3, 4.3,],
                                      [ 1.6, 3.4, 4.4,],
                                      [+80.7,  2.7, 3.7,],
                                      [+80.9,  2.9, 3.9,],
                                      [+80.1,  2.1, 3.1,],])            


    elif cell == "hex8":
        indexL = numpy.arange(16,20)
        indexN = numpy.arange(4,8)
        indexP = numpy.arange(12,16)
        n = 20
        m = 4
        DOF = 3

        A = numpy.array([[1.0,   0,   0,   0,],
                         [  0, 1.0,   0,   0,],
                         [  0,   0, 1.0,   0,],
                         [  0,   0,   0, 1.0,]])
        Cv = numpy.array([[ 0, +1, 0,],
                          [ 0, 0, +1,],
                          [ +1, 0, 0,],])
        Zv = numpy.zeros([3,3])
        C = numpy.vstack( (numpy.hstack((Cv, Zv, Zv, Zv)),
                           numpy.hstack((Zv, Cv, Zv, Zv)),
                           numpy.hstack((Zv, Zv, Cv, Zv)),
                           numpy.hstack((Zv, Zv, Zv, Cv)) ) )

        fieldT = numpy.array([[ 4.1, 2.1, 3.1,],
                              [ 4.2, 2.2, 3.2,],
                              [ 4.3, 2.3, 3.3,],
                              [ 4.4, 2.4, 3.4,],
                              [ 4.5, 2.5, 3.5,],
                              [ 4.6, 2.6, 3.6,],
                              [ 4.7, 2.7, 3.7,],
                              [ 4.8, 2.8, 3.8,],
                              [ 4.9, 2.9, 3.9,],
                              [ 4.0, 2.0, 3.0,],
                              [ 4.1, 2.1, 3.1,],
                              [ 4.2, 2.2, 3.2,],
                              [ 4.5, 2.5, 3.5,],
                              [ 4.6, 2.6, 3.6,],
                              [ 4.7, 2.7, 3.7,],
                              [ 4.8, 2.8, 3.8,],
                              [-4.4, 2.4, 3.4,],
                              [-4.6, 2.6, 3.6,],
                              [-4.8, 2.8, 3.8,],
                              [-4.0, 2.0, 3.0,],])

        if testCase == "stick":
            fieldTIncr = numpy.array([[0.1, 2.1, 1.1,],
                                      [0.2, 2.2, 1.2,],
                                      [0.3, 2.3, 1.3,],
                                      [0.4, 2.4, 1.4,],
                                      [0.5, 2.5, 1.5,],
                                      [0.6, 2.6, 1.6,],
                                      [0.7, 2.7, 1.7,],
                                      [0.8, 2.8, 1.8,],
                                      [0.9, 2.9, 1.9,],
                                      [0.0, 2.0, 1.0,],
                                      [1.1, 3.1, 2.1,],
                                      [1.2, 3.2, 2.2,],
                                      [0.5, 2.5, 1.5,],
                                      [0.6, 2.6, 1.6,],
                                      [0.7, 2.7, 1.7,],
                                      [0.8, 2.8, 1.8,],
                                      [-12.266666666667,  3.6, 4.6,],
                                      [-12.066666666667,  5.4, 2.4,],
                                      [-16.866666666667,  2.2, 8.2,],
                                      [-17.666666666667, 10.0, 2.0,],])
        elif testCase == "slip":
            fieldTIncr = numpy.array([[ 1.1, 2.1, 0.1,],
                                      [ 1.2, 2.2, 0.2,],
                                      [ 1.3, 2.3, 0.3,],
                                      [ 1.4, 2.4, 0.4,],
                                      [ 1.5, 2.5, 0.5,],
                                      [ 1.6, 2.6, 0.6,],
                                      [ 1.7, 2.7, 0.7,],
                                      [ 1.8, 2.8, 0.8,],
                                      [ 1.9, 2.9, 0.9,],
                                      [ 1.0, 2.0, 0.0,],
                                      [ 1.1, 2.1, 0.1,],
                                      [ 1.2, 2.2, 0.2,],
                                      [ 1.5, 2.9, 0.7,],
                                      [ 1.6, 2.8, 0.5,],
                                      [ 1.7, 2.9, 0.8,],
                                      [ 1.8, 2.7, 0.9,],
                                      [-1.4, 2.4, 0.4,],
                                      [-1.6, 2.6, 0.6,],
                                      [-1.8, 2.8, 0.8,],
                                      [-1.0, 2.0, 0.2,],])
        elif testCase == "open":
            fieldTIncr = numpy.array([[ 1.1, 2.1, 0.1,],
                                      [ 1.2, 2.2, 0.2,],
                                      [ 1.3, 2.3, 0.3,],
                                      [ 1.4, 2.4, 0.4,],
                                      [ 1.5, 2.5, 0.5,],
                                      [ 1.6, 2.6, 0.6,],
                                      [ 1.7, 2.7, 0.7,],
                                      [ 1.8, 2.8, 0.8,],
                                      [ 1.9, 2.9, 0.9,],
                                      [ 1.0, 2.0, 0.0,],
                                      [ 1.1, 2.1, 0.1,],
                                      [ 1.2, 2.2, 0.2,],
                                      [ 2.5, 1.9, 0.8,],
                                      [ 2.6, 1.8, 0.9,],
                                      [ 2.7, 1.7, 1.0,],
                                      [ 2.8, 1.6, 1.1,],
                                      [+20.4, 2.4, 0.4,],
                                      [+20.6, 2.6, 0.6,],
                                      [+20.8, 2.8, 0.8,],
                                      [+20.0, 2.0, 0.2,],])          

    # ------------------------------------------------------------------
    A /= lengthScale**2

    fieldTpdt = (fieldT + fieldTIncr)
    print "fieldTpdt",fieldTpdt
    fieldTpdtFault = globalToFault(fieldTpdt[indexL,:], C)

    tractionShearMag = (fieldTpdtFault[:,0]**2 + fieldTpdtFault[:,1]**2)**0.5
    tractionNormal = fieldTpdtFault[:,2]
    print "tractionShearMag",tractionShearMag
    print "tractionNormal",tractionNormal

    tractionRheologyFault = numpy.zeros((m, DOF))
    if testCase != "open":
        tractionRheologyShear = -0.6 * tractionNormal
        tractionRheologyNormal = tractionNormal

        tractionRheologyFault[:,0] = tractionRheologyShear * fieldTpdtFault[:,0]/tractionShearMag
        tractionRheologyFault[:,1] = tractionRheologyShear * fieldTpdtFault[:,1]/tractionShearMag
        tractionRheologyFault[:,2] = tractionRheologyNormal
    print "tractionRheologyFault",tractionRheologyFault

    tractionRheology = faultToGlobal(tractionRheologyFault, C)
    tractionInternal = fieldTpdt[indexL,:]
    tractionResidual = tractionRheology - tractionInternal
    print "tractionRheology",tractionRheology
    print "tractionInternal",tractionInternal
    print "tractionResidual",tractionResidual

    residual = numpy.zeros(fieldT.shape)
    residual[indexN,:] = +numpy.dot(A, tractionResidual)
    residual[indexP,:] = -numpy.dot(A, tractionResidual)
    residual[indexL,:] = numpy.dot(A, tractionInternal) * (fieldTpdt[indexP,:] - fieldTpdt[indexN,:])

    print "numpy.dot(A, tractionInternal)",numpy.dot(A, tractionInternal)
    print "fieldTpdt[indexP,:] - fieldTpdt[indexN,:]",fieldTpdt[indexP,:] - fieldTpdt[indexN,:]

    print "residual \n",printdata(residual)
