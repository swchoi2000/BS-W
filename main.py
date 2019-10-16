# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 15:10:55 2019

@author: NickolasT
"""
import os
from leafSolver import *
import translator as tr
import grapher as gr
import mapper as mp
import weightSolver as ws
import matplotlib.pyplot as plt
import numpy as np
import time

globalTime = time.time()
# This is the main file, excute this python script
# Include any necessary function calls here for testing purposes

# =====================PARAMETERS=====================
fileDelim = "\t"
nLvls = 3                           # Number of layers to use
xMin= -4.875                        # Minimum x-coord
xMax = 4.875                        # Maximum x-coord
xStep = 0.25                        # x descretization
searchStep = 1E-2                   # Step size for search values
weightSearchStep = 1E-2                # Threshold of weight convergence
maxWeight = 3                       # Maximum allowed weight
# ==============CONTROLS======================
DO_LEAF_POSITION_SOLVE = True       # True=Run solver.py methods to optimize leaf positions
                                    # (this only needs to be ran once)
DO_LEAF_EXPORT = True              # Use export functions from translator.py
                                    # to export leaf and weight info
DO_WEIGHT_SOLVE = True              # True=Run weightSolver.py method to optimize weight

rowEx = 4                           # Row index to display 1D fluence and leaf positions
tableWidth = 20                     # Char width of table entries
logFilename = "LOG.txt"
# Angles to check
angles = np.array(range(1,179))
#angles = ["160","320"]

global logFile
logFile = open(logFilename,'w')

if DO_LEAF_POSITION_SOLVE:
    tempTime = time.time()
    idealFluencesArray = []
    optimFluenceArray = []
    solverUpArray = []
    solverLowArray = []
    solverAvgArray = []
    for a,angle in enumerate(angles):
        angleTime = time.time()
        angle = str(angle)
        abs_file_path="C:/Users/staff/Documents/optimalFluence[2019-08-02]"
        filename = abs_file_path+"/eclipse178Field/Field "+angle+".optimal_fluence"     # File containing target fluence values
        expectedfilename = abs_file_path+"/ANM178Field/Field "+angle+".optimal_fluence"   
    
        # =========DO STUFF===================
        params = (nLvls, xMin, xMax, xStep,searchStep, maxWeight, weightSearchStep)
        # Read in the field data
        out = tr.importFluence(filename)
        optimalFluence = out[-1]
        xSize = out[1]
        ySize = out[2]
        
        # Simplify the field data by a factor of 4
        # Results in 10 rows = 10 leaf rows
        normFluence = tr.normalizedKernel(optimalFluence,4)
        nRows = len(normFluence)
        
        
        # Read in field data of expected solution
        out = tr.importFluence(expectedfilename)
        expectedFluence = out[-1]
        expectedFluence = tr.normalizedKernel(expectedFluence,4)
        
        # Instantiate a solver object
        mySolver1 = leafSolver(params)
        mySolver2 = leafSolver(params)
        mySolver3 = leafSolver(params)
        
        # Perform different solving algorithms 
        solve1, ideal, x, time1 = mySolver1.solveExtend2PeakAndTrough(normFluence,.05,True)
        solve2, ideal, x, time2 = mySolver2.solveExtend2PeakAndTrough(normFluence,.05,False)
        solve3, ideal, x, time3 = mySolver3.solveExtend2PeakWeightLimit(normFluence)
        
        idealFluencesArray.append(ideal)
        optimFluenceArray.append(expectedFluence)
        solverUpArray.append(mySolver1)
        solverLowArray.append(mySolver2)        
        solverAvgArray.append(mySolver3)
        
        # Calculate errors and graph for algorithm2
        fig, sumError1, sumErrorOver1,_ = gr.mapFluence(x,ideal,solve1,(0,1))        # Graph the fluence patterns
    #    fig.savefig(filename.replace(".","_")+"_"+str(nLvls)+"lvls.png")
        gr.mapFluence1D(x,ideal[rowEx],x,solve1[rowEx])      # Graph a row of the fluence
        gr.drawLeaves(mySolver1,rowEx)
        # Calculate for other
        _,sumError0, sumErrorOver0,_ = gr.diffMap(ideal,expectedFluence)
        _,sumError2, sumErrorOver2,_ = gr.diffMap(ideal,solve2)
        _,sumError3, sumErrorOver3,_ = gr.diffMap(ideal,solve3)
        
        # Print results
        print("===============================\nANGLE: %d degrees"% (int(angle)/len(angles)*360))
        print("===ALGORIGTHM===".ljust(tableWidth)
              +"===Weight===".ljust(tableWidth)
              +"===TimeElapsed===".ljust(tableWidth)
              +"===TotalError===".ljust(tableWidth)
              +"===TotalOver===".ljust(tableWidth))
        print("OPTIMAL".ljust(tableWidth)
              +("----").ljust(tableWidth)
              +("----").ljust(tableWidth)
              +("%.4E" % sumError0).ljust(tableWidth)
              +("%.4E" % sumErrorOver0).ljust(tableWidth))
        print("2PeakUpper".ljust(tableWidth)
              +("%.4f" % mySolver1.weight).ljust(tableWidth)
              +("%.5f" % time1 + "s").ljust(tableWidth)
              +("%.4E" % sumError1).ljust(tableWidth)
              +("%.4E" % sumErrorOver1).ljust(tableWidth))
        print("2PeakLower".ljust(tableWidth)
              +("%.4f" % mySolver2.weight).ljust(tableWidth)
              +("%.5f" % time2 + "s").ljust(tableWidth)
              +("%.4E" % sumError2).ljust(tableWidth)
              +("%.4E" % sumErrorOver2).ljust(tableWidth))
        print("2PeakAvg".ljust(tableWidth)
              +("%.4f" % mySolver3.weight).ljust(tableWidth)
              +("%.5f" % time3 + "s").ljust(tableWidth)
              +("%.4E" % sumError3).ljust(tableWidth)
              +("%.4E" % sumErrorOver3).ljust(tableWidth))
        t = (time.time()-angleTime,(time.time()-angleTime)*(len(angles)-a))
        print("AngleTime: %f s\tEstTime: %f s" % t)
        
        # ============WRITE TO LOG FILE=============
        logFile.write("===============================\nANGLE: %d degrees\n" % (int(angle)/len(angles)*360))
        logFile.write("===ALGORITHM===".ljust(tableWidth)
              +("----").ljust(tableWidth)
              +"===TimeElapsed===".ljust(tableWidth)
              +"===TotalError===".ljust(tableWidth)
              +"===TotalOver===".ljust(tableWidth)+"\n")
        logFile.write("OPTIMAL".ljust(tableWidth)
              +("----").ljust(tableWidth)
              +("----").ljust(tableWidth)
              +("%.4E" % sumError0).ljust(tableWidth)
              +("%.4E" % sumErrorOver0).ljust(tableWidth)+"\n")
        logFile.write("2PeakUpper".ljust(tableWidth)
              +("%.4f" % mySolver1.weight).ljust(tableWidth)
              +("%.5f" % time1 + "s").ljust(tableWidth)
              +("%.4E" % sumError1).ljust(tableWidth)
              +("%.4E" % sumErrorOver1).ljust(tableWidth)+"\n")
        logFile.write("2PeakLower".ljust(tableWidth)
              +("%.4f" % mySolver2.weight).ljust(tableWidth)
              +("%.5f" % time2 + "s").ljust(tableWidth)
              +("%.4E" % sumError2).ljust(tableWidth)
              +("%.4E" % sumErrorOver2).ljust(tableWidth)+"\n")
        logFile.write("2PeakAvg".ljust(tableWidth)
              +("%.4f" % mySolver3.weight).ljust(tableWidth)
              +("%.5f" % time3 + "s").ljust(tableWidth)
              +("%.4E" % sumError3).ljust(tableWidth)
              +("%.4E" % sumErrorOver3).ljust(tableWidth)+"\n")
        logFile.write("AngleTime: %f s\n" % t[0])
        # =============Export results=============
        if DO_LEAF_EXPORT:        
            fig.savefig("algo2_field_"+angle+"_"+str(nLvls)+"lvls.png")
    t = (time.time() - tempTime)
    print("*******LEAF SOLVE COMPLETED***********")
    print("Stage Time: %.4f s" % t)
    print("***************************")
    
    logFile.write("*******LEAF SOLVE COMPLETED***********\n")
    logFile.write("Stage Time: %.4f s\n" % t)
    logFile.write("***************************\n")
if DO_WEIGHT_SOLVE:
    tempTime = time.time()
    # Creates mask for fluence calculation to be of a cylinder (circle cross-section)
    circ = mp.createCircleMap(4.7,(xMin,xMax),(xMax,xMin), xSize, ySize)[0]
    # Import Tmr map array describing basic beam characteristics
    TmrMap = mp.importTmrMap("Tmr_Map.txt")[0]
    # Convert angles array from string to float
    paramAngles = []
    for i in range(len(angles)):
        paramAngles.append(float(angles[i])/len(angles)*360)
    # Tuple containing shared parameters of solvers for calculating 3D fluence volume
    solverParams = (xMin, xMax, xStep, ySize, xSize, circ, paramAngles,TmrMap)
    # Get the ideal volume0 and the expected experimental volume2
    
    # Calculate the expected results for comparison
    volume0 = mp.mapFluences(idealFluencesArray, solverParams)  # Ideal 3D fluence
    volume2 = mp.mapFluences(optimFluenceArray, solverParams)  # Optimal 3D fluence solution
    # Create weight array from weights found by solver objects
    weights = []
    for w in range(len(solverAvgArray)):
        weights.append(solverAvgArray[w].weight)
    # Do the weight optimization
    finalSolverArray,finalWeights = ws.genetic(volume0,[solverUpArray,solverLowArray,solverAvgArray],weights,solverParams,logFile)
    # Evaluate results of weight optimization
    fluencesArray = mp.getSolverArrayFluencesWithWeights(finalSolverArray,finalWeights)

    volume1 = mp.mapFluences(fluencesArray,solverParams)
    # Draw the final resulting volume
    diffMax = np.amax(np.abs(volume1-volume0))
    fluenceMax = np.amax(np.vstack([volume1,volume0]))
    for i in range(nRows):
        fig,_,_,_ = gr.mapFluenceREL(x,volume0[i],volume1[i],(0,np.amax(volume1[i])),diffMax,fluenceMax)
        fig.savefig("row"+str(i)+"_"+str(nLvls)+"lvls[FINAL].png")
    if DO_LEAF_EXPORT:    
        for w,mySolve in enumerate(finalSolverArray):
            filenameLeafOut = str(angles[w])+"_"+str(nLvls)+"lvls[FINAL].leaves"       # File to export leaf positions
            mySolve.weight = finalWeights[w]
            tr.outputLeaves(filenameLeafOut, mySolve,fileDelim)
    # Display the errors of the experimental
    print("Error: %f" % (mp.diffVolume(volume0,volume1)[1]))
    print("ExpectedError: %f" % (mp.diffVolume(volume0,volume2)[1]))
    logFile.write("Error: %f\n" % (mp.diffVolume(volume0,volume1)[1]))
    logFile.write("ExpectedError: %f\n" % (mp.diffVolume(volume0,volume2)[1]))
    
    print("*******WEIGHT SOLVE COMPLETED***********")
    print("Stage Time: %.4f s" % (time.time() - tempTime))
    print("***************************")
    logFile.write("*******WEIGHT SOLVE COMPLETED***********\n")
    logFile.write("Stage Time: %.4f s\n" % (time.time() - tempTime))
    logFile.write("***************************\n")
print("*******COMPLETED***********")
print("Total Time: %.4f s" % (time.time() - globalTime))
print("***************************")
logFile.write("*******COMPLETED***********\n")
logFile.write("Total Time: %.4f s\n" % (time.time() - globalTime))
logFile.write("***************************\n")
logFile.close()