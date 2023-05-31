#================================================================
# Import necessay libaries
#================================================================
import numpy as np
import scipy as sp
import scipy.optimize as spopt
import mkbs
import matplotlib.pyplot as plt
import latexutils as lutil


strVersion  = "Hedge Optimizer - 9.4.2014, python 3: 2023"
bTrace      = False
bTracePlus  = False

#================================================================
# This function tests whether s is a number - for processing the data
#================================================================

def is_number(s):
    try:
        temp = float(s)
    except:
        return (False, None)
    return (True, temp) 

def rosen(x):
    """The Rosenbrock function"""
    return sum(100.0*(x[1:]-x[:-1]**2.0)**2.0 + (1-x[:-1])**2.0)


class HedgeOptimizer:
    def __init__(self):

# ===================================================================
# Data Structure
# --------------
# MyAssetDir: Holds Put and other hedge Assets, both existing ones and new ones. For existing
#             ones the parameters are given, for new ones the Tool can put some of the parameters and 
#             leave others to be optimised
#             Nomeclatura "P" given put "P1", etc puts to be optimized (necissitates a tweak to mkbs routine)
#
# MyRequest = dict(): Here we store existing and new derivatives in the form
# MyRequest["1"] = [Instrument, Time, Strike]
# Entries as follows:
#       Instrument = "C","P","F","Z","B"
#       Time   Number or [Low, Max] in first case uses number in second case optimeses between low and high
#       Strike Number or [Low, Max] in first case uses number in second case optimeses between low and high
# MyRequest["Value"] = Number optimises for a derivative spend of Number in Total
# typical entry:
# -------------
#  val = psymObj[("P", Time, AssetPrice, Strike)]
#  there are four possible open parameters:
#  val [Number of instuments]: can be given, or free in form [min, max], with min or max also == None
#  Time: can be given, or free in form [min, max], with min or max also == None
#  AssetPrice <--- exogenically given via trading grid
#  Strike: can be given, or free in form [min, max], with min or max also == None
#
# MyParameterMap:       
# MyReverseParameterMap: Data which maps open Parameters to vector and back (Parameter Map to Vector and Reverse back to structure)
#
# MyTradingGridLia: This determines the Target Grid to be optimised including the relative weight:
#             example MyTradingGridLia[0.8,0.02,0.18,"pi"] = (2400, 0.08) where
#                     0.8: relative equity level (multiplicative)
#                     0.02: relative Interst Rate change (additive)
#                     0.02: relative Volatility level (additive)
#                     "pi": quantity to optimize (in {"pi","delta","gamma","rho","vega"})
#                     2400: liability value
#                     0.08: weigth for penalty 
#
# ===================================================================

        self.iCTracePenalty = -1
        self.strPenaly = ""
        self.MyAssetDir = dict()
        self.MyParameterMap = dict()
        self.MyReverseParameterMap = dict()
        self.MyRequest = dict()
        self.MyParameters = []
        self.MyLimits = []
        self.MyTradingGridLia = dict()
        self.Sx0 = None
        self.rZero= None
        self.SigmaZero = None
        self.bValueCondition = False
        self.fOutput = open("OutHO.txt","w")  
        self.LatexUtil = lutil.LatexUtils()

#================================================================
# The function ReturnColorDict - is used to have thw right colors for the bars
# pyplot is using RGB, hence the numbers (a,b,c) refer to the corresponding Red
# Green and Blue Values, which need to be between [0,1]
#================================================================

    def ReturnColorDict(self, AssetDir):
# MinMaxColor means that the colors go from the fist entry eg [0] to the last entry [1]
        MinMaxCol = { "B": ((0.8,0.8,0.8),(0.8,0.8,0.8)),
                      "F": ((1.0,1.0,0),(1.0,1.0,0)),
                      "P": ((0.,0.,0.9),(0.,0.9,0)),
                      "C": ((1.,0.,0.), (1., 1. ,0.))}

        iEntries = 0.
        MaxNr = dict()
        ActNr = dict()
        myColDir = dict()
        for i in AssetDir.keys():
            iEntries += 1.
            Letter = i[0][0]
            if Letter in MaxNr.keys():
                MaxNr[Letter] +=1
            else:
                MaxNr[Letter] =1
                ActNr[Letter] =0
        
        for i in AssetDir.keys():
            Letter = i[0][0]
            alpha = float(ActNr[Letter])/MaxNr[Letter]
            ActNr[Letter] +=1
            a = MinMaxCol[Letter][0]
            b = MinMaxCol[Letter][1]
            c = []
            for k in range(len(a)):
                c.append((1.-alpha) * a[k] + (alpha) * b[k])
            myColDir[i] = c
        return(myColDir, iEntries)
    
#================================================================
# Here we enter the starting point of the economy
#================================================================


    def InputStdEconomy(self,equitylevel, rate, vola):
        self.Sx0 = equitylevel
        self.rZero= rate
        self.SigmaZero = vola
        print("Std Economy set: S0 %10.2f, Rate %10.2f %%, Vola %.2f %%" % (self.Sx0, 100.*self.rZero,  100* self.SigmaZero))

        
#================================================================
# Here the trading grid is entered. For the moment being this is hard
# coded. An entry self.MyTradingGridLia[0.7,0.01,0.02,"pi"] = (9833.117, 1.)
# needs to be interpreted as follows
#   0.7: equity level (multiplicative)
#   0.01: interest rate shift applied eg + 1%
#   0.02: volatility shift applied eg +2%
#   9833.117: The resulting liability value "Value(k)"
#   1.: the weight to be applied for the optimisation "w_k"
# Note:
#   The optimiser used the folliwing penalty
#   \sum_k (F(k,x) -Value(k))^2 * w_k != min
#================================================================

    def InputTradingGrid(self):
        print("Simple Implementation of Trading Grid Input - eg hard coded")
        print("We also input key-economics here ....")
        print("Use 1.1.2012 - only partial input")
        print("Use 20.5%% foe vola and t=10yr irate at mid point")
        self.InputStdEconomy(1289.09,0.02151,0.205)
        self.MyTradingGridLia[0.7,0.0,0.0,"pi"] = (9833.117, 1.) 
        self.MyTradingGridLia[0.8,0.0,0.0,"pi"] = (8109.542, 1.) 
        self.MyTradingGridLia[0.9,0.0,0.0,"pi"] = (6714.820, 0.5) 
        self.MyTradingGridLia[1.0,0.0,0.0,"pi"] = (5671.411, 10000.) 
        self.MyTradingGridLia[1.1,0.0,0.0,"pi"] = (4945.008, 0.5) 
        self.MyTradingGridLia[1.2,0.0,0.0,"pi"] = (4476.287, 1.) 
        self.MyTradingGridLia[1.3,0.0,0.0,"pi"] = (4194.888, 1.) 
        self.MyTradingGridLia[1.0,0.0,0.0,"delta"] = (6.756*self.Sx0, 0) 
        self.MyTradingGridLia[1.0,0.0,0.0,"gamma"] = (0.005749*self.Sx0, 0.) 
        self.MyTradingGridLia[1.0,0.0,0.0,"rho"] = (28.119*1.e4, 0.0) 
        self.MyTradingGridLia[1.0,0.0,0.0,"vega"] = (280.237*1.e2, 0.0) 

# The following is only used to format the output

        strOut ="\n Trding Grid \n" 
        for i in self.MyTradingGridLia.keys():
            strOut +=" -- Key %s --> Value %s \n" % (str(i), str(self.MyTradingGridLia[i]))
        self.fOutput.write(strOut)
        self.LatexTradingGrid = ""
        AssetDir = self.MyTradingGridLia.keys()
        AssetDir=sorted(AssetDir)
        StrData = []
        NameEntriesX = ["$\\frac{S_t}{S_0}$","$\\Delta r$","$\\Delta \\nu$","Type","Value","$w_k$"]
        NameEntriesY = []
        iC = 0
        for i in AssetDir:
            iC += 1
            NameEntriesY.append(str(iC))
            X  = []
            for k in range(len(i)):
                # print i[k], type(i[k])
                # print "X:",X
                if type(i[k]) is list: 
                    X.append(str(i[k]))
                    # print "is list"
                else:
                    (a,b) = is_number(i[k])
                    # print a, b
                    if a:
                        X.append(self.LatexUtil.nicenumber(b,2)[0])
                    else:
                        X.append(str(i[k]))
            Values = self.MyTradingGridLia[i]
            X.append(self.LatexUtil.nicenumber(Values[0],2)[0])
            X.append(self.LatexUtil.nicenumber(Values[1],2)[0])
            StrData.append(X)
        matrix = np.matrix(StrData)
        Header ="\\section{Trading Grid} \n"
        self.LatexTradingGrid += self.LatexUtil.LatexPrintMatrix(matrix, Header, NameEntriesX, NameEntriesY,details = False, colors=False, nicenumbers=False)

#===============================================================
# Here the input to the optimisation task is defined
#================================================================


    def InputOptimisationTask(self,Instrument, Time, Strike, Amount):
        if bTrace:
            print("Input to opt Task - empty")
            print("Instrument", Instrument)
            print("Time", Time)
            print("Strike", Strike)
            print("Amount", Amount)
        if Instrument == "Value":
           Name = Instrument
        else:
            Name = Instrument + "%d" % (len(self.MyRequest))
        self.MyRequest[Name] = (Instrument, Time, Strike, Amount)
        if bTrace:
            print(self.MyRequest)

# ==================================================================
# In the following we have to map the free parameters of the optimisation
# task to the vector which the optimiser eats
# Concretely we need to define a map and the inicial conditions
# ==================================================================

    def PrepareOptimisationTask(self):
       if bTrace:
            print("Preparing Optimisation Task")
# ==================================================================
#
#   1. Prepare Map
#   2. Set Inicial Conditions
#
# ===================================================================
       self.MyParameterMap = dict()
       self.MyReverseParameterMap = dict()
       self.MyParameters = []
       self.MyLimits = []

       for i in self.MyRequest.keys():
           # print "\n Summary key", i

           # print "Parameter Map", self.MyParameterMap
           # print "Reverse Parameter Map", self.MyReverseParameterMap
           # print "Start Paraneters ", self.MyParameters
           # print "Limits", self.MyLimits 

           RequestedValues = self.MyRequest[i]
           if i == "Value":
               if bTrace:
                   print("++ Value set at ", RequestedValues[3])
               self.bValueCondition = RequestedValues[3]
           else:
               for k in range(len(RequestedValues)):
                   # print RequestedValues[k], type(RequestedValues[k])
                   if type(RequestedValues[k]) is list:
                       if bTrace:
                           print("Found Boundary Condition",i,"Entry ",k,"-->", RequestedValues[k])
                           # a) Set Start Value as new entry
                       (low, high) = RequestedValues[k]
                       if low == None:
                           startvalue = 0.8 * high
                       elif high == None:
                           startvalue = 1.2 * low
                       else:
                           startvalue = 0.5 * (high + low)
                       entrynr = len(self.MyParameters)
                       self.MyParameters.append(startvalue)
                       if bTrace:
                           print("Creating Entry nr %d with startvalue %f" % (entrynr, startvalue))
                           # b) Set boundaries for new entry
                       self.MyLimits.append((low, high))
                       if bTrace:
                           print(" ... with limits ", low, " and ", high)
                           # c) Update Parameter and Reverse Parameter Maps

                       self.MyParameterMap[(i,k)] = entrynr
                       self.MyReverseParameterMap[entrynr] = (i,k)
       if bTrace:
            print ("**** Ending PrepareOptimisationTask ***")
            print ("\n Summary")

            print ("Parameter Map", self.MyParameterMap)
            print ("Reverse Parameter Map", self.MyReverseParameterMap)
            print ("Start Paraneters ", self.MyParameters)
            print ("Limits", self.MyLimits)
            
            print ("try Creating Asset Dir ")
            ad = self.CreateLocalAssetDict(self.MyParameters)
            print (ad)

# ===================================================================
# For a given x vector from the optimised we need to construct the respective
# asset directory via the maps defined before
# Note that the coding is not particularly efficient - but has been kept for
# readability purposes...
# ===================================================================

    def CreateLocalAssetDict(self,x):
        locDict = dict()
        for i in self.MyRequest.keys():
            RequestedValues = self.MyRequest[i]
            if i != "Value":
                LocPar = []
                if i[0] == "P":
                    for k in range(len(RequestedValues)):
                   # print RequestedValues[k], type(RequestedValues[k])
                       if type(RequestedValues[k]) is list:
                           entrtyNr = self.MyParameterMap[(i,k)]
                           value = x[entrtyNr]
                           if bTracePlus:
                               print("Found Boundary Condition",i,"Entry ",k,"-->", RequestedValues[k], "replace by ", value)
                           LocPar.append(value)
                       else:
                           LocPar.append(RequestedValues[k])
                    mkbs.AddPut(locDict, LocPar[3] , LocPar[1], self.Sx0, LocPar[2]*self.Sx0)
                if i[0] == "C":
                    for k in range(len(RequestedValues)):
                   # print RequestedValues[k], type(RequestedValues[k])
                       if type(RequestedValues[k]) is list:
                           entrtyNr = self.MyParameterMap[(i,k)]
                           value = x[entrtyNr]
                           if bTracePlus:
                               print("Found Boundary Condition",i,"Entry ",k,"-->", RequestedValues[k], "replace by ", value)
                           LocPar.append(value)
                       else:
                           LocPar.append(RequestedValues[k])
                    mkbs.AddCall(locDict, LocPar[3] , LocPar[1], self.Sx0, LocPar[2]*self.Sx0)
                if i[0] == "F":
                    for k in range(len(RequestedValues)):
                   # print RequestedValues[k], type(RequestedValues[k])
                       if type(RequestedValues[k]) is list:
                           entrtyNr = self.MyParameterMap[(i,k)]
                           value = x[entrtyNr]
                           if bTracePlus:
                               print("Found Boundary Condition",i,"Entry ",k,"-->", RequestedValues[k], "replace by ", value)
                           LocPar.append(value)
                       else:
                           LocPar.append(RequestedValues[k])
                    mkbs.AddFuture(locDict, LocPar[3] , LocPar[1], self.Sx0, LocPar[2]*self.Sx0)
                if i[0] == "B":
                    for k in range(len(RequestedValues)):
                   # print RequestedValues[k], type(RequestedValues[k])
                       if type(RequestedValues[k]) is list:
                           entrtyNr = self.MyParameterMap[(i,k)]
                           value = x[entrtyNr]
                           if bTracePlus:
                               print("Found Boundary Condition",i,"Entry ",k,"-->", RequestedValues[k], "replace by ", value)
                           LocPar.append(value)
                       else:
                           LocPar.append(RequestedValues[k])
                    mkbs.AddCash(locDict, LocPar[3] , None, None, None)


        return(locDict)
# ===================================================================
# This defines the target  / penalty function for the optimiser
# There are two task involved - 1. To define the asset directory
#  and 2. to caluculate the penalty for the given vector x
#  Note there are two types of constraints 
#  a) The maximal spend for option see a) below and
#  b) The difference of values in a stressed environment and the
#     Greeks respectively see b) below
# ===================================================================


    def TargetFunction(self,x):
        if self.iCTracePenalty<0:
            self.strPenaly = ""
        if self.iCTracePenalty:
            self.strPenaly += "\n Iteration Penalty %d \n -----------------------\n" % (self.iCTracePenalty)
        AssetDict = self.CreateLocalAssetDict(x)
        penalty = 0.
        (fTotal0, delta0, gamma0, rho0, vega0) = mkbs.ValueOptions(AssetDict, 1.*self.Sx0, self.rZero, self.SigmaZero, symOutFile = None,bChatter = bTrace)
        target0, weight0 = self.MyTradingGridLia[(1., 0., 0., "pi")]
        if self.iCTracePenalty:
            self.strPenaly += "\n --- ValueO %f weight %f\n" % (target0, weight0)

        for (equitylevel, ilevel, volalevel, typestress) in self.MyTradingGridLia.keys():
            target, weight =  self.MyTradingGridLia[(equitylevel, ilevel, volalevel, typestress)]
            if self.iCTracePenalty:
                self.strPenaly += "\n --- Key %s - Value %f weight %f\n" % (str((equitylevel, ilevel, volalevel, typestress)),target, weight)

            # a) Calc Value of AssetDir
            if (equitylevel == 1.0) and (ilevel == 0) and (volalevel == 0) and typestress == "pi":
                 if self.bValueCondition:
                     # c) Do Value Condition on Top
                     penalty += weight * (self.bValueCondition-fTotal0)**2
                     if self.iCTracePenalty:
                         self.strPenaly += "\n ---- Value Penalty after %f - increment (weigth %f- FloatValue %f  TarValueO %f)\n" % (penalty, weight, self.bValueCondition, fTotal0)

            else:
                     (fTotalStress, deltaStress, gammaStress, rhoStress, vegaStress) = mkbs.ValueOptions(AssetDict, equitylevel*self.Sx0, self.rZero+ilevel, self.SigmaZero+volalevel, symOutFile = None,bChatter = bTrace)
            # b) Calc Penalty
                     if typestress == "pi":
                         penalty += weight * ((fTotalStress-fTotal0)-(target - target0))**2
                         if self.iCTracePenalty:
                             self.strPenaly += "\n ---- pi Penalty after %f - increment (weigth %f- FloatValue %f  TarValueO %f) piStress %f\n" % (penalty, weight, fTotalStress-fTotal0, target - target0,fTotalStress)

                     elif typestress == "delta":
                         penalty += weight * (deltaStress-target)**2
                         if self.iCTracePenalty:
                             self.strPenaly += "\n ---- delt Penalty after %f - increment (weigth %f- FloatValue %f  Target %f) \n" % (penalty, weight, target, deltaStress)                

                     elif typestress == "gamma":
                         penalty += weight * (gammaStress-target)**2
                     elif typestress == "rho":
                         penalty += weight * (rhoStress-target)**2
                     elif typestress == "vega":
                         penalty += weight * (vegaStress-target)**2
                         if self.iCTracePenalty:
                             self.strPenaly += "\n ---- delt Vega after %f - increment (weigth %f- FloatValue %f  Target %f) \n" % (penalty, weight, target, vegaStress)                

        if self.iCTracePenalty:
            self.iCTracePenalty -= 1
        return(penalty)
# ==================================================================
#
#   1. Read Parameters and update Asset Directory to Reflect Input
#   2. Calculate the Values for all Entries of the Trading Grid
#   3. Calulate Penalty
#
# ===================================================================
        

    # def ImputStartoint(self):
    #     if bTrace:
    #         print"Init Start"
    #     self.x0 = np.array([1.3, 0.7, 0.8, 1.9, 1.2])

    def Optimize(self):
        strOut = "Optimisation Tsk : \n"
        self.latexOut = "\\section{Optimizer Task} \n \\begin{verbatim} \n Requested Portfolio \n -------------------- \n"
        for i in self.MyRequest.keys():
            strOut +=" -- Key %s --> Value %s \n" % (str(i), str(self.MyRequest[i]))
            self.latexOut+=" -- Key %s --> Value %s \n" % (str(i), str(self.MyRequest[i]))
#   1. Read Parameters and update Asset Directory to Reflect Input
#   2. Calculate the Values for all Entries of the Trading Grid
        self.PrepareOptimisationTask()
        print(self.MyLimits)
        print("Dimension of problem %d" % (len(self.MyParameters)))
        # self.ImputStartoint()
        print("Do Optimize")
#   3. Calulate Penalty
        res = spopt.minimize(self.TargetFunction, self.MyParameters, method='L-BFGS-B', bounds =self.MyLimits, tol =1.e-12, options={'maxiter': 10000 , 'disp': True})
        # res = spopt.minimize(self.TargetFunction, self.MyParameters, method='TNC', bounds =self.MyLimits, tol =1.e-12, options={'maxiter': 10000 , 'disp': True})
        print(res)
        self.latexOut+=" \n \n Optimizer Result\n================== \n " + str(res) +"\n \\end{verbatim}\n"
        strOut += "\n\n Optimiser Output: \n %s \n Asset Dir: \n" % (str(res)) 
        self.AssetDict = self.CreateLocalAssetDict(res.x)
        for i in self.AssetDict.keys():
            strOut +=" -- Key %s --> Value %s \n" % (str(i), str(self.AssetDict[i]))
        strOut +="\n Value of Option Portfolio \n"
        self.fOutput.write(strOut)

# ==================================================================
#
#   The following two routines prepare the output including figures
#
# ===================================================================


        
    def Output(self):
        print("Doing Output")
        print("Asset Dir:")
        print(self.AssetDict)
        print("Valuation:")
        level = [ 1, 0.9, 0.8, 0.7]
        for i in level:
            print ("Valuation at level %f \n ------------------\n"%(i))
            (fTotal0, delta0, gamma0, rho0, vega0) = mkbs.ValueOptions(self.AssetDict, i*self.Sx0, self.rZero, self.SigmaZero, symOutFile = self.fOutput,bChatter = True)
        LevelSet = []
        LevelsTradingGridX = []
        LevelsTradingGridY = []
        LevelsTradingGridY2 = []

        for (equitylevel, ilevel, volalevel, typestress) in self.MyTradingGridLia.keys():
            if ilevel == 0 and volalevel == 0 and typestress == "pi":
                LevelSet.append(equitylevel)
        LevelSet=sorted(LevelSet)
        ilevel = 0
        volalevel = 0
        typestress = "pi"
        for equitylevel in LevelSet:
                LevelsTradingGridX.append(equitylevel)
                LevelsTradingGridY.append(self.MyTradingGridLia[(equitylevel, ilevel, volalevel, typestress)][0])
                (fTotal0, delta0, gamma0, rho0, vega0) = mkbs.ValueOptions(self.AssetDict, equitylevel*self.Sx0, self.rZero, self.SigmaZero, symOutFile = None,bChatter = False)
                LevelsTradingGridY2.append(fTotal0+ self.MyTradingGridLia[1.0,0.0,0.0,"pi"][0])
        print(LevelsTradingGridX)
        plt.figure(1)
        savename ="Fit.png"
        plt.plot(LevelsTradingGridX,LevelsTradingGridY,"-r",LevelsTradingGridX,LevelsTradingGridY2,"pc",markersize=12)
        plt.grid(True)
        plt.savefig(savename, dpi=300, orientation='landscape', papertype='a4')            
        plt.figure(2)
        for i in range(len(LevelsTradingGridX)):
            if LevelsTradingGridX[i] == 1.:
                k = i
        for i in range(len(LevelsTradingGridX)):

            LevelsTradingGridY2[i] = LevelsTradingGridY2[i] + LevelsTradingGridY[k] - LevelsTradingGridY2[k]
        savename ="FitShift.png"
        plt.plot(LevelsTradingGridX,LevelsTradingGridY,"-r",LevelsTradingGridX,LevelsTradingGridY2,"pc",markersize=12)
        plt.grid(True)
        plt.savefig(savename, dpi=300, orientation='landscape', papertype='a4')            

    def OutputLatex(self, psymFile=None, fileName = None, FullOutput = False, Pics = True):
        print("Doing Latex Output")
        iNrPic = 0
        locDigits = 2
        start = ""
        end   = ""
        if FullOutput:
            start = self.LatexUtil.LatexTokens["docmatter"][0]
            end = self.LatexUtil.LatexTokens["docmatter"][1]
        strOutput = start + self.LatexTradingGrid +"\n" + self.latexOut
        if psymFile == None:
            if fileName:
                tempname = fileName
            else:
                tempname = "testlatex"
            psymFile = open (tempname+".tex","w")
        print("Asset Dir:")
        AssetDir = self.AssetDict.keys()
        AssetDir=sorted(AssetDir)
        StrData = []
        NameEntriesX = ["Type","$t$","$S_t$","Strike","Number"]
        NameEntriesY = []
        iC = 0
        for i in AssetDir:
            iC += 1
            NameEntriesY.append(str(iC))
            X  = []
            print("i:",i)
            for k in range(len(i)):
                # print i[k], type(i[k])
                # print "X:",X
                if type(i[k]) is list: 
                    X.append(str(i[k]))
                    # print "is list"
                else:
                    (a,b) = is_number(i[k])
                    # print a, b
                    if a:
                        X.append(self.LatexUtil.nicenumber(b,locDigits)[0])
                    else:
                        X.append(str(i[k]))
            X.append(self.LatexUtil.nicenumber(self.AssetDict[i],locDigits)[0])
            StrData.append(X)
        print(StrData)
        matrix = np.matrix(StrData)
        Header ="\\section{Asset Directory} \n"
        strOutput += self.LatexUtil.LatexPrintMatrix(matrix, Header, NameEntriesX, NameEntriesY,details = False, colors=False, nicenumbers=False)

        Header ="\\section{Valuation} \n"
        StrData = []
        NameEntriesX = ["$\\pi$","$\\delta$","$\\gamma$","$\\rho$","$\\nu$"]
        NameEntriesY = []
        level = [ 1.3, 1.2, 1.1, 1, 0.9, 0.8, 0.7]
        for i in level:
            NameEntriesY.append("%6.2f"%(i))
            print("Valuation at level %f \n ------------------\n"%(i))
            # (fTotal0, delta0, gamma0, rho0, vega0) = mkbs.ValueOptions(self.AssetDict, i*self.Sx0, self.rZero, self.SigmaZero, symOutFile = self.fOutput,bChatter = True)
            X = mkbs.ValueOptions(self.AssetDict, i*self.Sx0, self.rZero, self.SigmaZero, symOutFile = self.fOutput,bChatter = True)
            StrData.append(X)
        matrix = np.matrix(StrData)
        strOutput += self.LatexUtil.LatexPrintMatrix(matrix, Header, NameEntriesX, NameEntriesY,details = False, colors=False, nicenumbers=True, digits = locDigits)

        for i in level:
            StrData = []
            Header ="\\subsection{Valuation Details at level %4.2f} \n" % (i)
            StrData = []
            NameEntriesX = ["Type","$t$","$S_t$","Strike","$\\pi$","$\\delta$","$\\rho$"]
            NameEntriesY = []
            XX = mkbs.ValueOptions(self.AssetDict, i*self.Sx0, self.rZero, self.SigmaZero, symOutFile = None,bChatter = False, longOutPut=True)

            iC = 0
            for i in AssetDir:
                iC += 1
                NameEntriesY.append(str(iC))
                X  = []
                print("i:",i)
                for k in range(len(i)):
                    # print i[k], type(i[k])
                    # print "X:",X
                    if type(i[k]) is list: 
                        X.append(str(i[k]))
                        # print "is list"
                    else:
                        (a,b) = is_number(i[k])
                        # print a, b
                        if a:
                            X.append(self.LatexUtil.nicenumber(b,locDigits)[0])
                        else:
                            X.append(str(i[k]))
                Values = XX[-1][i]
                for k in Values:
                    X.append(self.LatexUtil.nicenumber(k,locDigits)[0])
                StrData.append(X)
            print(StrData)
            matrix = np.matrix(StrData)
            strOutput += self.LatexUtil.LatexPrintMatrix(matrix, Header, NameEntriesX, NameEntriesY,details = False, colors=False, nicenumbers=False)
            




        LevelSet = []
        LevelsTradingGridX = []
        LevelsTradingGridY = []
        LevelsTradingGridY2 = []
        LevelsTradingGridY3 = []
        for (equitylevel, ilevel, volalevel, typestress) in self.MyTradingGridLia.keys():
            if ilevel == 0 and volalevel == 0 and typestress == "pi":
                LevelSet.append(equitylevel)
        LevelSet=sorted(LevelSet)
        ilevel = 0
        volalevel = 0
        typestress = "pi"
        latexpicout = ""
        for equitylevel in LevelSet:
                LevelsTradingGridX.append(equitylevel)
                LevelsTradingGridY.append(self.MyTradingGridLia[(equitylevel, ilevel, volalevel, typestress)][0])
                (fTotal0, delta0, gamma0, rho0, vega0) = mkbs.ValueOptions(self.AssetDict, equitylevel*self.Sx0, self.rZero, self.SigmaZero, symOutFile = None,bChatter = False)
                LevelsTradingGridY2.append(fTotal0+ self.MyTradingGridLia[1.0,0.0,0.0,"pi"][0])
                LevelsTradingGridY3.append(fTotal0)
        print(LevelsTradingGridX)

        plt.figure(1)
        savename =tempname+("%d.png" %(iNrPic))
        iNrPic +=1
        caption = "Adjusted Comparision"
        latexpicout +=  self.LatexUtil.LatexTokens["figure"][0]
        latexpicout +=  self.LatexUtil.LatexTokens["graphic"][0] %(savename)
        latexpicout +=  self.LatexUtil.LatexTokens["figure"][1] % (caption)
        plt.plot(LevelsTradingGridX,LevelsTradingGridY,"-r",LevelsTradingGridX,LevelsTradingGridY2,"pc",markersize=12)
        plt.grid(True)
        plt.savefig(savename, dpi=300, orientation='landscape', papertype='a4')            

        plt.figure(2)
        for i in range(len(LevelsTradingGridX)):
            if LevelsTradingGridX[i] == 1.:
                k = i
        for i in range(len(LevelsTradingGridX)):
            LevelsTradingGridY2[i] = LevelsTradingGridY2[i] + LevelsTradingGridY[k] - LevelsTradingGridY2[k]
        savename =tempname+("%d.png" %(iNrPic))
        iNrPic +=1
        caption = "Adjusted Comparision"
        latexpicout +=  self.LatexUtil.LatexTokens["figure"][0]
        latexpicout +=  self.LatexUtil.LatexTokens["graphic"][0] %(savename)
        latexpicout +=  self.LatexUtil.LatexTokens["figure"][1] % (caption)
        plt.plot(LevelsTradingGridX,LevelsTradingGridY,"-r",LevelsTradingGridX,LevelsTradingGridY2,"pc",markersize=12)
        plt.grid(True)
        plt.savefig(savename, dpi=300, orientation='landscape', papertype='a4')            
        
        locColDict,NrBars =  self.ReturnColorDict(self.AssetDict)
                      
        xBar = []
        hBar = []
        bBar = []
        colBar = []

        for equitylevel in LevelSet:
            h = 0.
            xOffset = -0.09 / 2.
            for k in AssetDir:
                
                locAD = dict()
                locAD[k] = self.AssetDict[k] 
                (fTotal0, delta0, gamma0, rho0, vega0) = mkbs.ValueOptions(locAD, equitylevel*self.Sx0, self.rZero, self.SigmaZero, symOutFile = None,bChatter = False)
                xBar.append(equitylevel + xOffset) 
                xOffset  += 0.09 / NrBars
                if fTotal0 >= 0:
                    hBar.append(fTotal0) 
                    bBar.append(h)
                    h += fTotal0
                else:
                    hBar.append(-fTotal0) 
                    bBar.append(h+fTotal0)
                    h += fTotal0
                colBar.append(locColDict[k])
                
        plt.figure(3) 
        # plt.title("Investment Grade Bonds")
        plt.bar(xBar, hBar, bottom = bBar, color =colBar, width = 0.09 / NrBars)
        #plt.hold(True)
        plt.plot(LevelsTradingGridX,LevelsTradingGridY3,"-",color = (0.8,0.8,0.8))
        plt.grid(True)

        savename =tempname+("%d.png" %(iNrPic))
        iNrPic +=1
        caption = "Decomposition P\&L"
        latexpicout +=  self.LatexUtil.LatexTokens["figure"][0]
        latexpicout +=  self.LatexUtil.LatexTokens["graphic"][0] %(savename)
        latexpicout +=  self.LatexUtil.LatexTokens["figure"][1] % (caption)
        plt.savefig(savename, dpi=300, orientation='landscape', papertype='a4')

        tMin = 100000.
        locAD = dict()
        for k in self.AssetDict.keys():
            t = k[1]
            if t !=None and k[0] != 'B':
                tMin = min(t, tMin)
        print("tMin = ", tMin)
        for k in self.AssetDict.keys():
            kTilde = []
            for i in range(len(k)):
                 kTilde.append(k[i])
            if k[1] != None:
                kTilde[1] = k[1]+0.0005-tMin
            locAD[tuple(kTilde)] = self.AssetDict[k]
        
        xVect = []
        yVect = []
        for i in level:
            XX = mkbs.ValueOptions(locAD, i*self.Sx0, self.rZero, self.SigmaZero, symOutFile = None,bChatter = False, longOutPut=True)
            xVect.append(i)
            yVect.append(XX[0])

        plt.figure(4) 
        plt.title("Response Function at t=%.2f" % (tMin))
        plt.plot(xVect,yVect)
        plt.grid(True)
        savename =tempname+("%d.png" %(iNrPic))
        iNrPic +=1
        plt.savefig(savename, dpi=300, orientation='landscape', papertype='a4')


        if Pics:
            strOutput += latexpicout
        psymFile.write(strOutput + "\n" + end)

        
def main():
    task = 1
    print("Starting %s" % (strVersion))
    myopt = HedgeOptimizer()
    myopt.InputTradingGrid()
    myopt.InputOptimisationTask("Value", None, None, 250.)
    # myopt.InputOptimisationTask("P", 1, [0.6, 0.7],[0,1000])
    # myopt.InputOptimisationTask("P", 1, [0.7,0.8],[0,1000])
    if task == 2: 
        myopt.InputOptimisationTask("P", [1, 2] , [0.7,0.9],[0,1000])
        myopt.InputOptimisationTask("P", [0.5, 2] , [0.5,0.6],[0,200])
        myopt.InputOptimisationTask("C", [0.8, 1.5], [1.15,1.4],[-2,2])
        myopt.InputOptimisationTask("F", None, 1. ,[-500,500])
        myopt.InputOptimisationTask("B", None, None ,[0,250])
    else:
        myopt.InputOptimisationTask("P", [1, 2] , [0.5,0.6],[0,200])
        myopt.InputOptimisationTask("P", [4, 5] , [0.5,0.6],[0,200])
        myopt.InputOptimisationTask("P", [9, 10] , [0.5,0.6],[0,200])
        myopt.InputOptimisationTask("P", [1, 2] , [0.8,0.9],[0,200])
        myopt.InputOptimisationTask("P", [4, 5] , [0.8,0.9],[0,200])
        myopt.InputOptimisationTask("P", [9, 10] , [0.8,0.9],[0,200])
        myopt.InputOptimisationTask("C", [1, 1.2], [1.2,1.4],[-2,2])
        myopt.InputOptimisationTask("C", [4, 4.2], [1.2,1.4],[-2,2])
        myopt.InputOptimisationTask("C", [9, 9.2], [1.2,1.4],[-2,2])
    #   myopt.InputOptimisationTask("F", None, 1. ,[-500,500])
        myopt.InputOptimisationTask("B", None, None ,[0,250])
    #    myopt.InputOptimisationTask("P", 0.25, [0.6, 0.7],[0,1000])
    #    myopt.InputOptimisationTask("P", 0.25, [0.7,0.8],[0,1000])
    #    myopt.InputOptimisationTask("P", 0.25, [0.8,0.9],[0,1000])
    #    myopt.InputOptimisationTask("P", 0.25, [1.1,1.2],[-200,0])


    myopt.Optimize()
    myopt.Output()

    myopt.OutputLatex(fileName = "task4", FullOutput = True, Pics = True)

    print("Penalty Calculation:")
    print(myopt.strPenaly)

if __name__=="__main__":
    main()
