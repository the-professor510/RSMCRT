import numpy as np
import plotInverseClass
import matplotlib.pyplot as plt

folderName = "data/inverse/"

#fileNames = ["inverseAttempt1,1.dat", "inverseAttempt1,2.dat", "inverseAttempt1,3.dat", "inverseAttempt1,4.dat", "inverseAttempt1,5.dat", "inverseAttempt1,6.dat"]
fileNames = ["inverse.dat"]
musList = []
muaList = []
hggList = []
nList = []
bestGuessIndxList = []
bestMusList = []
bestMuaList = []
besthggList = []
bestnList = []
bestErrorList = []

for i in range(len(fileNames)):
    filename = folderName + fileNames[i]

    readInverse = plotInverseClass.plotInverseClass()
    mus, mua, hgg, n, error, bestGuessIndx, bestMus, bestMua, besthgg, bestn, bestError = readInverse.read_1D_Detector(filename)
    
    musList.append(mus)
    muaList.append(mua)
    hggList.append(hgg)
    nList.append(n)
    bestGuessIndxList.append(bestGuessIndx)
    bestMusList.append(bestMus)
    bestMuaList.append(bestMua)
    besthggList.append(besthgg)
    bestnList.append(bestn)
    bestErrorList.append(bestError)

    """
    fig, ax = plt.subplots(2,2, figsize=(8,8), layout="constrained")

    ax[0][0].scatter(mus, error)
    ax[0][0].scatter(mus[bestGuessIndx], error[bestGuessIndx], color = "red", label = "Best Guess")
    ax[0][0].set_xlabel("$\\mu_{s}$ $[distance^{-1}]$")
    ax[0][0].set_ylabel("Error [Arb. Units]")
    ax[0][0].legend()

    ax[0][1].scatter(mua, error)
    ax[0][1].scatter(mua[bestGuessIndx], error[bestGuessIndx], color = "red", label = "Best Guess")
    ax[0][1].set_xlabel("$\\mu_{a}$ $[distance^{-1}]$")
    ax[0][1].set_ylabel("Error [Arb. Units]")
    ax[0][1].legend()

    ax[1][0].scatter(hgg, error)
    ax[1][0].scatter(hgg[bestGuessIndx], error[bestGuessIndx], color = "red", label = "Best Guess")
    ax[1][0].set_xlabel("$g$ [Unitless]")
    ax[1][0].set_ylabel("Error [Arb. Units]")
    ax[1][0].legend()

    ax[1][1].scatter(mus*(1-hgg), error)
    ax[1][1].scatter((mus*(1-hgg))[bestGuessIndx], error[bestGuessIndx], color = "red", label = "Best Guess")
    ax[1][1].set_xlabel("$\\mu_{s}' = \\mu_{s}(1-g)$ $[distance^{-1}]$")
    ax[1][1].set_ylabel("Error [Arb. Units]")
    ax[1][1].legend()

    plt.show()
    #"""

    readInverse.plot1D(mus, error, bestGuessIndx, "mus")
    readInverse.plot1D(mua, error, bestGuessIndx, "mua")
    readInverse.plot1D(np.sqrt(3*mua*(mua+mus*(1-hgg))), error, bestGuessIndx, "mueff diffusion approximation")
    #readInverse.plot1D(hgg, error, bestGuessIndx, "hgg")
    #readInverse.plot1D(n, error, bestGuessIndx, "n")
    #readInverse.plot1D(mus*(1-hgg), error, bestGuessIndx, "mus'")

    readInverse.plot2D(mus*(1-hgg), mua, error, bestGuessIndx, xName = "mus'", yName="mua")
    #readInverse.plot2D(mus, mua, error, bestGuessIndx, xName = "mus", yName="mua")
    #readInverse.plot2D(mus, hgg, error, bestGuessIndx, xName = "mus", yName="hgg")
    #readInverse.plot2D(mus, n, error, bestGuessIndx, xName = "mus", yName="n")

    #readInverse.plot2D(mua, hgg, error, bestGuessIndx, xName = "mua", yName="hgg")
    #readInverse.plot2D(mua, n, error, bestGuessIndx, xName = "mua", yName="n")

    #readInverse.plot2D(hgg, n, error, bestGuessIndx, xName = "hgg", yName="n")

"""
musList = np.array(musList)
muaList = np.array(muaList)
hggList = np.array(hggList)
nList = np.array(nList)
bestGuessIndxList = np.array(bestGuessIndxList)
bestMusList = np.array(bestMusList)
bestMuaList = np.array(bestMuaList)
besthggList = np.array(besthggList)
bestnList = np.array(bestnList)
bestErrorList = np.array(bestErrorList)

fig, ax = plt.subplots(2,2, figsize=(8,8), layout="constrained")
edges = np.arange(0.5,7.5,1)
minList = [0,0,0,0,4.5,4.8]
maxList = [50,50,50,15,5.5,5.2]
#ax[0][0].set_yscale('log')
ax[0][0].stairs(bestMusList, edges, label = "Best Guess", linewidth = 1.5)
ax[0][0].stairs(minList,edges, label = "Maximum Bound", linewidth = 1.5)
ax[0][0].stairs(maxList,edges, label = "Minimum Bound", linewidth = 1.5)
ax[0][0].hlines(5.0,0.5,6.5, label = "True Value", color = "tab:red", linestyle = "--")
ax[0][0].set_xlim(0.5,6.5)
#ax[0][0].set_ylim(0,50.0)
ax[0][0].set_ylabel("$\\mu_{s}$ $[distance^{-1}]$")
ax[0][0].set_xlabel("iteration")
ax[0][0].legend(loc="upper right")


minList = [0,0,1.75,1.9,1.925,1.95]
maxList = [50,10,2.25,2.1,2.075,2.03]
#ax[0][1].set_yscale('log')
ax[0][1].stairs(bestMuaList, edges, label = "Best Guess", linewidth = 1.5)
ax[0][1].stairs(minList,edges, label = "Maximum Bound", linewidth = 1.5)
ax[0][1].stairs(maxList,edges, label = "Minimum Bound", linewidth = 1.5)
ax[0][1].hlines(2.0,0.5,6.5, label = "True Value", color = "tab:red", linestyle = "--")
ax[0][1].set_xlim(0.5,6.5)
#ax[0][1].set_ylim(0,50.0)
ax[0][1].set_ylabel("$\\mu_{a}$ $[distance^{-1}]$")
ax[0][1].set_xlabel("iteration")
ax[0][1].legend(loc="upper right")

minList = [0,0,0,0,0.3,0.35]
maxList = [1,1,1,1,0.5,0.46]
ax[1][0].stairs(besthggList, edges, label = "Best Guess", linewidth = 1.5)
ax[1][0].stairs(minList,edges, label = "Maximum Bound", linewidth = 1.5)
ax[1][0].stairs(maxList,edges, label = "Minimum Bound", linewidth = 1.5)
ax[1][0].hlines(0.4,0.5,6.5, label = "True Value", color = "tab:red", linestyle = "--")
ax[1][0].set_xlim(0.5,6.5)
#ax[1][0].set_ylim(0,1.0)
ax[1][0].set_ylabel("$g$ [Unitless]")
ax[1][0].set_xlabel("iteration")
ax[1][0].legend(loc="upper right")

minList = [0,0,1.5,2.5,2.8,2.9]
maxList = [50,10,4.5,3.5,3.1,3.1]
#ax[1][1].set_yscale('log')
ax[1][1].stairs(bestMusList*(1-besthggList),edges, label = "Best Guess", linewidth = 1.5)
ax[1][1].stairs(minList,edges, label = "Maximum Bound", linewidth = 1.5)
ax[1][1].stairs(maxList,edges, label = "Minimum Bound", linewidth = 1.5)
ax[1][1].hlines(3.0,0.5,6.5, label = "True Value", color = "tab:red", linestyle = "--")
ax[1][1].set_xlim(0.5,6.5)
#ax[1][1].set_ylim(0,50.0)
ax[1][1].set_ylabel("$\\mu_{s}' = \\mu_{s}(1-g)$ $[distance^{-1}]$")
ax[1][1].set_xlabel("iteration")
ax[1][1].legend(loc="upper right")

plt.show()
#"""
