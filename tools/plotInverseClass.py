import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import sys

class plotInverseClass:
       
    def plot1D(self, x, error, bestIndex, xName = "optical Prop"):
             
        
        fig = plt.figure(1)
        ax1 = fig.add_subplot()

        ax1.scatter(x, error, c = error)
        ax1.scatter(x[bestIndex], error[bestIndex], color = "red", label = "Best Guess")
        ax1.set_xlabel(xName)
        ax1.set_ylabel("Error (Arb. Units)")
        ax1.legend()
        plt.show()
        
    def plot2D(self, x, y, error, bestIndex, xName = "optical Prop A", yName = "optical Prop B"):
        fig = plt.figure(1)
        #ax1 = fig.add_subplot(projection='3d')
        ax1 = fig.add_subplot()

        
        #ax1.scatter(x,y, error, c=error)
        ax1.scatter(x,y, c=error)
        ax1.scatter(x[bestIndex], y[bestIndex], color = "blue", label = "Best Guess", marker = "x")
        ax1.set_xlabel(xName)
        ax1.set_ylabel(yName)
        
        
        # Plot the postior distribution and some samples
        fig, ax = plt.subplots(subplot_kw={"projection": "3d", "computed_zorder": False})
        try:
            ax.plot_trisurf(x,y, error, antialiased=True)
            ax.scatter(x[bestIndex], y[bestIndex], error[bestIndex], color = "blue", label = "Best Guess", marker = "x")
            ax.set_xlabel(xName)
            ax.set_ylabel(yName)
        except:
            print("error: couldn't perform surface plot")
            
        plt.show()
        
    #def plot3D(self, x, y, z, error, bestIndex):
        


    def read_1D_Detector(self, filename):
        data = np.fromfile(file=filename, dtype=np.float64, sep="")
        
        
        bestGuessIndx, bestMus, bestMua, bestHgg, bestn, bestError = data[:6]
        bestGuessIndx = int(bestGuessIndx)
        
        mus, mua, hgg, n, error = data[6:].reshape(-1, 5).T
            
        return mus, mua, hgg, n, error, bestGuessIndx, bestMus, bestMua, bestHgg, bestn, bestError

                    
           
