import ROOT
import pandas as pd

def ntuple_to_pd(path,runs,channels):
    mylist=[]
    for run in runs:
        aux=[]
        df1 = ROOT.RDataFrame("ntuple",path+"run%i_NTuple.root" %run)
        Qdf1 = ROOT.RDataFrame("charge",path+"run%i_NTuple.root" %run)
        diccionario = df1.AsNumpy()
        Qdiccionario=Qdf1.AsNumpy()
        aspd = pd.DataFrame(diccionario,columns=None)
        Qaspd=pd.DataFrame(Qdiccionario,columns=None)
        all_pd=pd.DataFrame.join(aspd,Qaspd)
        
        for ch in channels:
            aux.append(all_pd[all_pd["ch"]==ch])
        mylist.append(aux)
    
    return mylist