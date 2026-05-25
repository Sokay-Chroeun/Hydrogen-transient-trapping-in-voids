import assembly
import part
import visualization
import operator
import odbAccess
import time
import mesh #pour lancer en noGUI, modif 8 mars 2013
import os
import numpy as np
import subprocess
from abaqus import *
from abaqusConstants import *
from caeModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()

def produit(UneListe,DeuxListe):
  compteur = 0
  accumulateur = []
  if len(UneListe) == len(DeuxListe):
    largeur = len(UneListe)
    while compteur < largeur:
      accumulateur.append(UneListe[compteur]*DeuxListe[compteur])
      compteur += 1
  return accumulateur

# somme de tous les termes d'une liste
def somme(UneListe):
  compteur = 0
  accumulateur = 0
  largeur = len(UneListe)
  while compteur < largeur:
    accumulateur += UneListe[compteur]
    compteur += 1
  return accumulateur

# moyenne ponderee de la liste echantillon
def stat_moyenne( echantillon, poids=[]) :
  if poids==[]:
    moyenne = somme( echantillon ) / len(echantillon)
  else:
    moyenne = somme( produit(echantillon, poids) ) / somme(poids)
  return moyenne

# variance ponderee de la liste echantillon
def stat_variance( echantillon, poids=[]) :
  variance = stat_moyenne( produit(echantillon,echantillon), poids) -  stat_moyenne( echantillon, poids)**2
  return variance

# variance ponderee de la liste echantillon
def stat_ecart_type( echantillon , poids=[]) :
  ecart_type = sqrt( stat_variance( echantillon, poids ) )
  return ecart_type

def postTrait(jobName):
    odbName=jobName
    monPath = os.getcwd()+'/'
    os.chdir(monPath)
    openOdb = session.openOdb(name=odbName+'.odb',readOnly=True)
    stepName = session.odbs[odbName+'.odb'].steps.keys()
    currentOdb = session.odbs[odbName+'.odb']
    odbStep = currentOdb.steps
    Stat_filename = odbName+'-Qtot.txt'
    f = open(Stat_filename, 'w')
    #f.write('openOdb: %r \n' % (odbName))
    f.write('time CL CT CH\n' )
    for CurrentStepName in stepName:
      myOdb = odbStep[CurrentStepName]
      frame = myOdb.frames
      tempsInitialStep = myOdb.totalTime -1. #step-0 : mise en place
      for currentFrame in frame:
        currentFieldOutputs = currentFrame.fieldOutputs
        vol = currentFieldOutputs['IVOL'].values #volume de ts les pts de Gauss
        CL = currentFieldOutputs['SDV1'] # CL
        CT = currentFieldOutputs['SDV2'] # CT2
        CLValues = CL.values # CL
        CTValues = CT.values # CT
        GPVol=[] # vol d'un pt de Gauss
        GPCL=[] # CL sur  pt de Gauss
        GPCT=[] # CT2 sur pt de Gauss
        for i in range(0,len(vol)): # sauvegarde des donnees
          CurrentVol = vol[i]
          GPVol.append(CurrentVol.data)
          CurrentCL = CLValues[i]
          GPCL.append(CurrentCL.data)
          CurrentCT = CTValues[i]
          GPCT.append(CurrentCT.data)
        QL = somme(produit(GPCL, GPVol))
        QT = somme(produit(GPCT, GPVol))
        VolumeGP = somme(GPVol)
        strOuput = str(tempsInitialStep+currentFrame.frameValue+1)+" "+str(QL)+" "+str(QT)+" "+str(QT+QL)+"\n"
        f.write(strOuput)
    f.close()
    print('post processing finished')
    return(0)

postTrait('Job-1')