import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from math import exp, log as ln
from tabulate import tabulate

# constants
q = 1.602 * 10**(-19)
k = 1.381 * 10**(-23)
T = 298.15

def calcVoc(Isc, Io):
    Voc = ((k*T)/q) * ln((Isc/Io) + 1)
    return Voc  

def calcCurrent(Isc, Io, Vd, Rp):
    I = Isc - Io*(exp(38.9*Vd) -1) - (Vd/Rp)
    return I

def calcV(I, Vd, Rs):
    V = Vd - (I * Rs);
    return V

def calcP(V, I):
    return V*I

def getIArray(VdArray, Isc, Io, Rp):
    IArray = []
    for Vd in VdArray:
        I = calcCurrent(Isc, Io, Vd, Rp)
        IArray.append(round(I, 4))
    return IArray

def getVshArray(VdArr, IArr, Rp, Rs, numCells):
    VshArray = []
    # print(((numCells -1)/numCells) * V)
    for(V, I) in zip(VdArr, IArr):
        Vsh = (((numCells - 1)/numCells)*(V*numCells)) - (I * (Rp + Rs))
        # if Vsh > 0:
        VshArray.append(round(Vsh, 4))
        
    return VshArray

class PVModule:
    def __init__(self, Isc, Io, Rp, Rs, numCells):
        self.numCells = numCells
        self.Isc = Isc
        self.Io = Io
        self.Rp = Rp
        self.Rs = Rs
        self.Voc = calcVoc(Isc, Io)
        self.VdArr = [x * 0.001 for x in range(0, round(self.Voc*1000))]
        self.IArr = getIArray(self.VdArr, self.Isc, self.Io, self.Rp)
        self.VshArr = getVshArray(self.VdArr, self.IArr, self.Rp, self.Rs, self.numCells)

    def plot(self):
        fig, ax1 = plt.subplots()
        ax1.set_xlim(0, 68)
        ax1.set_ylim(0, 6.3)
        ax1.set_xlabel("Voltage (V)")
        ax1.set_ylabel("Current (A)")
        color = 'tab:red'
        ax1.plot(self.VshArr, self.IArr, '--', label="Nth Cell Shaded", color=color)
        ax1.yaxis.set_major_locator(ticker.MaxNLocator(20))
        ax1.xaxis.set_major_locator(ticker.MaxNLocator(20))
        color = 'tab:blue'
        ax1.plot(np.multiply(self.VdArr, self.numCells), self.IArr, '-', label="Unshaded", color=color)
        plt.legend(loc="lower left")
        plt.title("Unshaded & Nth Cell Shaded IV Curves")
        plt.show()


class PVCell:
    def __init__(self, cell, Isc, Io, Rp, Rs):
        self.cell = cell
        self.Isc = Isc
        self.Io = Io
        self.Rp = Rp
        self.Rs = Rs
        self.Voc = calcVoc(Isc, Io)
        self.VdArr = [x * 0.001 for x in range(0, round(self.Voc*1000))]
        self.IArr = getIArray(self.VdArr, self.Isc, self.Io, self.Rp)
        self.PArr = np.multiply(self.IArr, self.VdArr)
        self.Vmpp = self.VdArr[np.argmax(self.PArr)] # Max P Vmpp
        self.Impp = self.IArr[np.argmax(self.PArr)] # i am so clever
        self.Pmpp = self.Vmpp * self.Impp
        self.FF = self.Pmpp/(self.Voc * self.Isc)

    def printVoc(self):
        print("Voc = " + str(round(self.Voc, 4)))
    
    def printVdArr(self):
        print(self.VdArr)
    
    def printIArr(self):
        print(self.IArr)

    def printMPP(self):
        print("Vmpp = " + str(round(self.Vmpp, 4)))
        print("Impp = " + str(round(self.Impp, 4)))

    def printPmpp(self):
        print("Pmpp = " + str(round(self.Pmpp, 4)))
 
    def printFF(self):
        FF = self.Pmpp/(self.Voc * self.Isc)
        print("FF = " + str(round(FF, 4)))

    def printInfo(self):
        print("For cell " + self.cell)
        self.printVoc()
        self.printMPP()
        self.printPmpp()
        self.printFF()
    
    def plot(self):
        fig, ax1 = plt.subplots()
        color = 'tab:blue'
        ax1.set_xlabel("Voltage (V)")
        ax1.set_ylabel("Current (A)")
        ax1.plot(self.VdArr, self.IArr, '-', label="IV Curve", color=color)
        ax1.plot(self.Vmpp, self.Impp, 'ro', label="MPP")
        ax1.yaxis.set_major_locator(ticker.MaxNLocator(20))
        ax1.xaxis.set_major_locator(ticker.MaxNLocator(10))
        plt.legend(loc="center left")
        
        ax2 = ax1.twinx()
        color = 'tab:red'
        ax2.set_ylabel("Power (W)")
        ax2.plot(self.VdArr, self.PArr, '--', label="Power", color=color)
        ax2.yaxis.set_major_locator(ticker.MaxNLocator(10))

        fig.tight_layout()
   
        plt.legend(loc="lower left")
        plt.title("Cell "+ self.cell)
        plt.show()

if __name__ == "__main__":
    cellA = PVCell("A", 5, 6*(10**-11), 10, 0.001)
    cellB = PVCell("B", 6, 5*(10**-11), 12, 0.0025)
    cellC = PVCell("C", 5.5, 7*(10**-11), 8, 0.0015)
    module = PVModule(6, 5*(10**-11), 12, 0.0025, 102)


    data = [[cellA.cell, cellA.Voc, cellA.Vmpp, cellA.Impp, cellA.Pmpp, cellA.FF],
            [cellB.cell, cellB.Voc, cellB.Vmpp, cellB.Impp, cellB.Pmpp, cellB.FF],
            [cellC.cell, cellC.Voc, cellC.Vmpp, cellC.Impp, cellC.Pmpp, cellC.FF]]

    print(tabulate(data, headers=["PV Cell", "Voc (V)", "Vmpp (V)", "Impp (A)", "Pmpp (W)", "FF"]))

    cellA.plot()
    
    cellB.plot()

    cellC.plot()

    module.plot()
