import scipy as sp
import numpy as np
import math

class LatexUtils:
    def __init__(self):
        self.RelativePrecision= 4
        self.Digits = 4
        self.Eps = 1e-5            # Epsilon
        self.bChatter = False
        self.strRedMacro = "\\textred{%s}"
        self.strGreenMacro = "\\textgreen{%s}"
        self.LatexTokens ={"tabular"  :("\\begin{center} \n %s \\par \n \\begin{tabular}{%s} \n","\\end{tabular} \n \\end{center}"),
                           "smalltab" :("\\begin{center} \n %s \\par \n \\begin{scriptsize} \n \\begin{tabular}{%s} \n","\\end{tabular} \n  \\end{scriptsize} \n \\end{center}"),
                           "figure"   :("\\begin{figure} \n","\\caption{%s} \n \\end{figure}"),
                           "graphic"  :("\\includegraphics[width=0.9\\textwidth]{%s} \n",""),
                           "pagewidth":(105., "mm"),
                           "docmatter":("""\\documentclass[12pt]{article} \n \\usepackage[pdftex]{graphicx} \n \\setlength{\\paperwidth}{210mm} \n 
                                         \\setlength{\\paperheight}{297mm}%% \n \\usepackage{a4}%% \n \\textwidth=16cm%% \n \\textheight=23cm%% \n 
                                         \\oddsidemargin=0.0cm%% \n \\evensidemargin=0.0cm%% \n \\parindent=0mm \n \\begin{document} \n""","\\end{document}")        }                            

    def LatexPrintMatrix(self, matrix, Header, NameEntriesX, NameEntriesY, colors=False, nicenumbers=False, details=True, digits = None):
        strOut = "%% Automatically generated Matrix Output \n"
        (N,M) = matrix.shape
        locDigits = self.Digits
        if digits != None:
            locDigits = digits
        pLen = self.LatexTokens["pagewidth"][0] / (M+1.)
        FormatToken ="p{%.1f %s}" %(pLen, self.LatexTokens["pagewidth"][1])
        Format = ""
        for k in range(M+1):
            Format += FormatToken 
            if k != M:
                Format += "|"
        strOut += self.LatexTokens["smalltab"][0] % (Header,Format)
        Header = "%10s"  % ("")
        for k in range(M):
               Header += "& %8s " % (NameEntriesX[k])
        strOut += Header + "\\\\ \hline \n"
        for i in range(N):
            Line = "%s " % (NameEntriesY[i])
            for k in range(M):
                number = matrix[i,k]
                roundednumber = number
                if nicenumbers:
                    strnumber, roundednumber = self.nicenumber(number,locDigits)
                else:
                    strnumber = str(number)
                if colors:
                    if roundednumber<0:
                        strnumber = self.strRedMacro % (strnumber)
                    elif roundednumber>0:
                        strnumber = self.strGreenMacro % (strnumber)
                Line += "& " + strnumber 
            strOut += Line + "\\\\ \n"
        strOut += self.LatexTokens["smalltab"][1]
        if details:
           (eigenval, eigenvec) = np.linalg.eig(matrix)
           b= np.zeros((len(eigenval),len(eigenval)+1))
           b[:,0] = eigenval.T
           b[:,1:] = eigenvec
           NamesXX = []
           NamesYY = []
           NamesXX.append("Eigenval")
           for i in range(len(eigenvec)):
               NamesXX.append("Dim %d" % (i))
               NamesYY.append("EV %d" % (i))
           strOut += self.LatexPrintMatrix(b, "EV", NamesXX, NamesYY, colors=True, nicenumbers=True, details=False, digits=2)
        return strOut
    

    def nicenumber(self, number, digits):
# This small routine provides nicely formated numbers
         if np.isnan(number):
             return("NaN",0)
         if np.isinf(number):
             return("$\\pm \\infty$",0)
         roundednumber = number
         if number >= 0:
             strNumber=""
             sign =+1.
         else:
             strNumber="-"
             sign =-1.
             number = -number
         if self.RelativePrecision < 10 and number > 0:
             mypower = int(math.log(number) / math.log(10))
             number = round(number / 10.0**mypower, self.RelativePrecision) * 10.0**mypower + self.Eps
             roundednumber = sign*number
         if abs(roundednumber) < 10.0**(-digits):
             roundednumber =0
             strNumber=""             
         myexp = range(9,-1,-3)
         biggerNrExists = 0;
         fraction = number - int(number)
         for i in myexp:
             (partnumber, number) = divmod(number, 10**i)
             if i == 9 and partnumber > 0:
                 strNumber += str(int(partnumber))+"'"
                 biggerNrExists = 1;
             elif i == 0:
                 tempstr =str(int(partnumber))
                 if biggerNrExists == 1:
                     tempstr=tempstr.rjust(3)
                     tempstr=tempstr.replace(' ','0')
                 strNumber += tempstr
             else:
                 tempstr =str(int(partnumber))
                 if biggerNrExists == 1:
                     tempstr=tempstr.rjust(3)
                     tempstr=tempstr.replace(' ','0')
                     strNumber += tempstr+"'"
                 else:
                     if partnumber > 0:
                         biggerNrExists = 1
                         strNumber += tempstr+"'"
         if digits > 0:
             strNumber += '.'
             tempstr =str(int(10**digits * fraction))
             tempstr=tempstr.rjust(digits)
             tempstr=tempstr.replace(' ','0')
             strNumber += tempstr
         return(strNumber,roundednumber)

    
def main():
    data = [[-0.2234444, 0.673674563],[1234444, -4673647563]]
    matrix = sp.matrix(data)
    Header = "Example"
    NameEntriesX = ["X1", "X2"]
    NameEntriesY = ["Bla Y 1","Y2"]
  
    MyObj = LatexUtils()
    
    Types = []
    Types.append((False, False))
    Types.append((False, True))
    Types.append((True, False))
    Types.append((False, True))
    iC = 0

    for (a,b) in Types:
        print("Run %d \n \n"% (iC))
        print(MyObj.LatexPrintMatrix(matrix, Header, NameEntriesX, NameEntriesY, colors=a, nicenumbers=b))
        iC +=1

    strmatrixdata= [["abcdef","hhrt"],["abcde","s"]]
    matrix = sp.matrix(strmatrixdata)
    print(MyObj.LatexPrintMatrix(matrix, Header, NameEntriesX, NameEntriesY,details = False, colors=False, nicenumbers=False))



if __name__ == "__main__":
     main()
