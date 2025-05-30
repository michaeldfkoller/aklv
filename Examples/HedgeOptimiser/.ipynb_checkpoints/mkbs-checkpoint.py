# Use
# Prepare structure
# MyAssetDir = dict()
# Adding puts and zeros as follows
#         if a1 == "Z":
#              AddZero(MyAssetDir, p[iCount], a2, a3, a4)
#         else:
#              AddPut(MyAssetDir, p[iCount], a2, a3/1000., a4/1000.)
# Valuation
# ValueOptions(MyAssetDir, 0.8, 0.025, 0.25,symOutFile = symOut)

#================================================================
#  1. Importing cumulative density function for normal distribution
#     Name normcdfgen
#================================================================

import math
try:
# Try whether we have scipy library to use the function - otherwise do it yourself 
    import scipy.stats as sps
    normcdfgen =sps.norm.cdf
except:
    print("use own normdist cdf")
    def erfcc(x):
        """Complementary error function."""
        z = abs(x)
        t = 1. / (1. + 0.5*z)
        r = t * math.exp(-z*z-1.26551223+t*(1.00002368+t*(.37409196+
          t*(.09678418+t*(-.18628806+t*(.27886807+
          t*(-1.13520398+t*(1.48851587+t*(-.82215223+
          t*.17087277)))))))))
        if (x >= 0.):
            return r
        else:
            return 2. - r

    def normcdfgen(x, mu=0, sigma=1):
        t = x-mu;
        y = 0.5*erfcc(-t/(sigma*math.sqrt(2.0)));
        if y>1.0:
            y = 1.0;
        return y

    x = [0., 2.57]
    for i in x:
        print ("%.1f  --> %.1f" % (i, normcdfgen(i)))

import pickle

#================================================================
# 2. Here we define a bunch of functions to add various assets, namely
#     Put ("P"), Call("C"), Zeros("Z"), Future("F") and Cash ("B" [bank])
#     The structure of the dictonary to set up is always the same
#     val = psymObj[("P", Time, AssetPrice, Strike)]
#     where
#          "P": sybol of the instrument
#           val: the number of instruments
#================================================================


def AddPut(psymObj, Amount, Time, AssetPrice, Strike):
    val = 0.
    if ("P", Time, AssetPrice, Strike) in psymObj.keys():
        val = psymObj[("P", Time, AssetPrice, Strike)]
    psymObj[("P", Time, AssetPrice, Strike)] = val + Amount

def AddCall(psymObj, Amount, Time, AssetPrice, Strike):
    val = 0.
    if ("C", Time, AssetPrice, Strike) in psymObj.keys():
        val = psymObj[("C", Time, AssetPrice, Strike)]
    psymObj[("C", Time, AssetPrice, Strike)] = val + Amount


def AddZero(psymObj, Amount, Time, AssetPrice, Strike):
    val = 0.
    if ("Z", Time, AssetPrice, Strike) in psymObj.keys():
        val = psymObj[("Z", Time, AssetPrice, Strike)]
    psymObj[("Z", Time, AssetPrice, Strike)] = val + Amount

def AddFuture(psymObj, Amount, Time, AssetPrice, Strike):
    val = 0.
    if ("F", Time, AssetPrice, Strike) in psymObj.keys():
        val = psymObj[("F", Time, AssetPrice, Strike)]
    psymObj[("F", Time, AssetPrice, Strike)] = val + Amount

def AddCash(psymObj, Amount, Time, AssetPrice, Strike):
    val = 0.
    Time = 0
    AssetPrice = None
    Strike = None
    if ("B", Time, AssetPrice, Strike) in psymObj.keys():
        val = psymObj[("B", Time, AssetPrice, Strike)]
    psymObj[("B", Time, AssetPrice, Strike)] = val + Amount


def normcdf(x):
    return(normcdfgen(x))

#================================================================
# 3. In the following we have used b/s code for value options
#    including the respective Greeks
#================================================================

def BlackScholesEuro(Price, CallPut, AssetP, Strike, RiskFree, DividendY, Time, Volatility):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  BlackScholesEuro(Price, CallPut, AssetP, Strike, RiskFree, DividendY, Time, Volatility):
#% Computes the Black-Scholes European Call/Put Option Values based
#% on the following inputs:
#% CallPut           =       Call = 1, Put = 0
#% AssetP            =       Underlying Asset Price
#% Strike            =       Strike Price of Option
#% RiskFree          =       Risk Free rate of interest
#% Time              =       Time to Maturity
#% Volatility        =       Volatility of the Underlying
#% Please note that the use of this code is not restricted in anyway.
#% However, referencing the author of the code would be appreciated.
#% To run this program, simply use the function defined in the 1st line.
#% http://www.global-derivatives.com
#% info@global-derivatives.com
#% Kevin Cheng (Nov 2003)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    dt = Volatility * (Time) ** 0.5
    df = RiskFree - DividendY + 0.5 * Volatility ** 2           # Computes the drift term
    d1 = (math.log(AssetP / Strike) + df * Time) / dt                #Calculates the d1 term used in Black-Scholes
    d2 = d1 - dt                                                # Calculates the d2 term used in Black-Scholes
    alpha = math.exp(-DividendY * Time)
# The cumulative normal distribution functions for use in computing calls
    nd1 = normcdf(d1)
    nd2 = normcdf(d2)
# The cumulative normal distribution functions for use in computing puts
    nnd1 = normcdf(-d1)
    nnd2 = normcdf(-d2)

    if CallPut == 1:
    # Computes call price
        Price = AssetP * alpha * nd1 - Strike * math.exp(-RiskFree * Time) * nd2
    else:
  # Computes put price
        Price = Strike * math.exp(-RiskFree * Time) * nnd2 - AssetP * alpha * nnd1
    return([Price, CallPut, AssetP, Strike, RiskFree, DividendY, Time, Volatility])


#function [] = BlackScholesMertonEuro(AssetP, Strike, RiskFree, DividendY, Time, Volatility)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#% Computes the Black-Scholes-Merton European Call/Put Option Values based
#% on the following inputs:
#% AssetP            =       Underlying Asset Price
#% Strike            =       Strike Price of Option
#% RiskFree          =       Risk Free rate of interest
#% DividendY         =       Dividend Yield of Underlying
#% Time              =       Time to Maturity
#% Volatility        =       Volatility of the Underlying
#% Please note that the use of this code is not restricted in anyway.
#% However, referencing the author of the code would be appreciated.
#% To run this program, simply use the function defined in the 1st line.
#% http://www.global-derivatives.com
#% info@global-derivatives.com
#% Kevin Cheng (Nov 2003)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

#$dt = Volatility * \sqrt{Time}$;
#$df = RiskFree - DividendY + 0.5 * Volatility ^ 2$;            % Computes the drift term
#$d1 = (\log( AssetP / Strike) + df * Time ) / dt$;             % Calculates the d1 term used in Black-Scholes
#$d2 = d1 - dt$;                                                % Calculates the d2 term used in Black-Scholes

#% The cumulative normal distribution functions for use in computing calls
#nd1 = normcdf(d1);
#nd2 = normcdf(d2);
#% The cumulative normal distribution functions for use in computing puts
#nnd1 = normcdf(-d1);
#nnd2 = normcdf(-d2);

#% Computes call price
#CallPrice = AssetP * Exp(-DividendY * Time) * nd1 - Strike * Exp(-RiskFree * Time) * nd2
#% Computes put price
#PutPrice = Strike * Exp(-RiskFree * Time) * nnd2 - AssetP * Exp(-DividendY * Time) * nnd1#


def Greeks(Delta, Gamma, Theta, Rho1, Rho2, Vega, CallPut, AssetP, Strike, RiskFree, DividendY, Time, Volatility):
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#% Computes the Black-Scholes-Merton European Call/Put Greeks based
#% on the following inputs:
#% CallPut           =       For a call, the input is "1". Puts, use "0".
#% AssetP            =       Underlying Asset Price
#% Strike            =       Strike Price of Option
#% RiskFree          =       Risk Free rate of interest
#% DividendY         =       Dividend Yield of Underlying
#% Time              =       Time to Maturity
#% Volatility        =       Volatility of the Underlying
#% Please note that the use of this code is not restricted in anyway.
#% However, referencing to the global derivatives website would be appreciated.
#% To run this program, simply use the function defined in the 1st line.
#% http://www.global-derivatives.com
#% info@global-derivatives.com
#% Kevin Cheng (Nov 2003)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Pi = 3.141592653589
    dt = Volatility * (Time) ** 0.5
    df = RiskFree - DividendY + 0.5 * Volatility ** 2                       # Computes the drift term
    d1 = (math.log(AssetP / Strike) + df * Time) / dt                # Calculates the d1 term used in Black-Scholes
    d2 = d1 - dt                                                # Calculates the d2 term used in Black-Scholes

# The cumulative normal distribution functions for use in computing calls
    nd1 = normcdf(d1)
    nd2 = normcdf(d2)
# The cumulative normal distribution functions for use in computing puts
    nnd1 = normcdf(-d1)
    nnd2 = normcdf(-d2)
    nn1 = (1 / (2 * Pi) ** 0.5 * math.exp(-0.5 * d1 ** 2))

# Computes call greeks
    if CallPut == 1:
        Delta = nd1 * math.exp(-DividendY * Time)
        Gamma = (nn1 * math.exp(-DividendY * Time)) / (AssetP * dt)
        Theta = -((AssetP * nn1 * math.exp(-DividendY * Time) * Volatility) / 2 * (Time) ** 0.5) + (-DividendY * AssetP * nd1 * math.exp(-DividendY * Time)) - (RiskFree * Strike * math.exp(-RiskFree * Time) * nd2)
        Vega = AssetP * (Time) ** 0.5 * nn1 * math.exp(-DividendY * Time)
        Rho1 = Strike * Time * math.exp(-RiskFree * Time) * nd2
        Rho2 = -AssetP * math.exp(-DividendY * Time) * Time * nd1
    else:
# Computes put greeks
        Delta = (nd1 - 1) * math.exp(-DividendY * Time)
        Gamma = (nn1 * math.exp(-DividendY * Time)) / (AssetP * dt)
        Theta = -((AssetP * nn1 * Volatility * math.exp(-DividendY * Time)) / 2 * (Time) ** 0.5) - (DividendY * AssetP * nnd1 * math.exp(-DividendY * Time)) + (RiskFree * Strike * math.exp(-RiskFree * Time) * nnd2)
        Vega = AssetP * (Time) ** 0.5 * nn1 * math.exp(-DividendY * Time)
        Rho1 = -Strike * Time * math.exp(-RiskFree * Time) * nnd2
        Rho2 = AssetP * math.exp(-DividendY * Time) * Time * nnd1

    return([Delta, Gamma, Theta, Rho1, Rho2, Vega, CallPut, AssetP, Strike, RiskFree, DividendY, Time, Volatility])

#================================================================
# 4. The function ValueOptions() is the main valuation function 
#    with the input as per below. It eats an asset directory (psym0), 
#    of the form generated above and values it
#     Input:
#     St: Equity level at time of valuation
#     iRF: (risk free) valuation interest rate
#     fSigma: volatility
#     symOutFile: if not None - output is written on the respective filepointer
#     bChatter: if False the subroutie shuts up (good if used many times)
#     longOutPut: if True the values of the individual instruments is retured as a diectory
#
#================================================================

def ValueOptions(psymO, St, iRF, fSigma, symOutFile = None,bChatter = True, longOutPut = False):
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  BlackScholesEuro(Price, CallPut, AssetP, Strike, RiskFree, DividendY, Time, Volatility):
#% Computes the Black-Scholes European Call/Put Option Values based
#% on the following inputs:
#% CallPut           =       Call = 1, Put = 0
#% AssetP            =       Underlying Asset Price
#% Strike            =       Strike Price of Option
#% RiskFree          =       Risk Free rate of interest
#% Time              =       Time to Maturity
#% Volatility        =       Volatility of the Underlying
    if bChatter:
        print("Summary Options ")
        print("----------------")
    strText = ["Funds Level   %6.1f", \
                   "Risk Free     %6.3f", \
                   "Sigma         %6.3f"]
    fValues = [ St, iRF, fSigma ]
    if symOutFile:
        strTemp = "\n"
        symOutFile.write(strTemp+"\n")

    for i in range(len(strText)):
         if bChatter:
            print(strText[i] % (fValues[i]))
         if symOutFile:
             strTemp = strText[i] % (fValues[i])
             symOutFile.write(strTemp+"\n")
    if bChatter:
        print( "               ")
    ValueDict = dict()
    fValVect = 50 * [0]
    delta = 0.
    gamma = 0.
    rho = 0.
    vega = 0.
    fTotal = 0.
    for k in psymO.keys():
#% Here we go through all assets which we received and value them includin their Greeks  
      # print k, symO[k]
        u = k[0]
        DividendY = 0.
        if u == 'P':
            iType = 0
        else:
            iType = 1
        if k[1] != 0 and u != 'Z' and u != 'F' and u != "B":    
            res = psymO[k] * BlackScholesEuro(0, iType, St, k[3], iRF, DividendY, k[1], fSigma)[0]
            iTime = max(0,min(49,int(k[1] + 0.5)))
            fValVect[iTime] += res
            fTotal += res

        # Calc Greeks       Greeks(Delta, Gamma, Theta, Rho1, Rho2, Vega, CallPut, AssetP, Strike, RiskFree, DividendY, Time, Volatility)
            gr = Greeks(0., 0., 0., 0., 0., 0., iType, St, k[3], iRF, DividendY, k[1], fSigma)
            locdelta = psymO[k] * gr[0] * St
            delta += locdelta
            gamma += psymO[k] * gr[1] * St
            locrho = psymO[k] * gr[3]
            rho   += locrho 
            vega  += psymO[k] * gr[5]
        if u == 'Z':
            iTime = max(0,min(49,int(k[1] + 0.5)))             
            locdelta = 0
            locrho = - psymO[k] * iTime
            rho   -= psymO[k] * iTime
            res = psymO[k] * (1+iRF)**(-iTime)
            fValVect[iTime] += res
            fTotal += res 
        if u == 'F':
            # val = psymObj[("F", Time, AssetPrice, Strike)]
            res = psymO[k] * (St-k[3])
            # print psymO[k],St,k[3]
            iTime = 0
            fValVect[iTime] += res
            fTotal += res
            locdelta = psymO[k] * St
            locrho = 0
            delta += locdelta
        if u == 'B':
            # val = psymObj[("F", Time, AssetPrice, Strike)]
            res = psymO[k] 
            # print psymO[k],St,k[3]
            iTime = 0
            locdelta = 0
            locrho = 0
            fValVect[iTime] += res
            fTotal += res
        ValueDict[k] = (res, locdelta, locrho)

    strTemp =  " delta    = %10.1f \n" % (delta)
    strTemp += " gamma    = %10.1f \n" % (gamma)
    strTemp += " rho      = %10.1f \n" % (rho)
    strTemp += " vega     = %10.1f \n" % (vega)    
    strTemp += " total    = %10.1f \n" % (fTotal)

    strTemp += "\n  allocation over time: \n"
    for iTime in range(50):
        if fValVect[iTime] !=0:
            strTemp += " t %-2d:    = %10.1f \n" % (iTime,fValVect[iTime])

    strTemp += "\n Value by Instrument \n"
    for k in ValueDict.keys():
        strTemp += "%-40s  = %10.1f %10.1f %10.1f\n" % (str(k), ValueDict[k][0], ValueDict[k][1], ValueDict[k][2])

    if bChatter:
        print(strTemp)
    if symOutFile:
        symOutFile.write(strTemp+"\n")
    if longOutPut:
        return(fTotal, delta, gamma, rho, vega, ValueDict)
    return(fTotal, delta, gamma, rho, vega)

#================================================================
# 5. This main() is only used for testing puropses
# 
# to run type: "python mkbs.py"
#================================================================



def main():
    print( "test mkbs")
    myAssetDir = dict()
    AddPut(myAssetDir, 1000., 1., 100000., 80000.)
    testfile = open("test.txt","w")
    ValueOptions(myAssetDir, 100000., 0.025, 0.18, symOutFile = testfile, bChatter = True)
    ValueOptions(myAssetDir, 99000., 0.025, 0.18, symOutFile = testfile, bChatter = True)
    ValueOptions(myAssetDir, 101000., 0.025, 0.18, symOutFile = testfile, bChatter = True)
    ValueOptions(myAssetDir, 100000., 0.0251, 0.18, symOutFile = testfile, bChatter = True)   
    ValueOptions(myAssetDir, 100000., 0.0249, 0.18, symOutFile = testfile, bChatter = True)   
    ValueOptions(myAssetDir, 100000., 0.025, 0.19, symOutFile = testfile, bChatter = True)
    ValueOptions(myAssetDir, 100000., 0.025, 0.17, symOutFile = testfile, bChatter = True)
    myAssetDir = dict()
    AddFuture(myAssetDir, 1000., 1., 1., 1.)
    ValueOptions(myAssetDir, 1., 0.025, 0.18, symOutFile = testfile, bChatter = True)
    ValueOptions(myAssetDir, 0.8, 0.025, 0.18, symOutFile = testfile, bChatter = True)
    testfile.close()


if __name__=="__main__":
    main()
