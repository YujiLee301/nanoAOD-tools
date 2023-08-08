#! /usr/bin/env python3
# Sources:
#   https://cms-nanoaod-integration.web.cern.ch/integration/master-106X/mc106Xul18_doc.html#GenPart
#   https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/python/genparticles_cff.py
#   https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/plugins/LHETablesProducer.cc
from __future__ import print_function # for python3 compatibility
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.tools import hasbit


def getprodchain(part,genparts=None,event=None):
  """Print production chain recursively."""
  chain = "%3s"%(part.pdgId)
  imoth = part.genPartIdxMother
  while imoth>=0:
    if genparts:
      moth = genparts[imoth]
      chain = "%3s –> "%(moth.pdgId)+chain
      imoth = moth.genPartIdxMother
    elif event:
      chain = "%3s –> "%(event.GenPart_pdgId[imoth])+chain
      imoth = event.GenPart_genPartIdxMother[imoth]
  return chain
  

def getdecaychain(part,genparts,indent=0,depth=999):
  """Print decay chain recursively."""
  chain   = "%3s"%(part.pdgId)
  imoth   = part._index
  ndaus   = 0
  indent_ = len(chain)+indent
  for idau in range(imoth+1,len(genparts)): 
    dau = genparts[idau]
    if dau.genPartIdxMother==imoth: # found daughter
      if ndaus>=1:
        chain += '\n'+' '*indent_
      if depth>0:
        chain += " –> "+getdecaychain(dau,genparts,indent=indent_+4,depth=depth-1)
      else: # stop recursion
        chain += " –> %3s"%(dau.pdgId)
      ndaus += 1
  return chain
  

# DUMPER MODULE
class LHEDumper(Module):
  
  def __init__(self):
    self.nleptonic = 0
    self.ntauonic  = 0
    self.nevents   = 0
  
  def analyze(self,event):
    """Dump gen information for each gen particle in given event."""
    print("\n%s Event %s %s"%('–'*10,event.event,'–'*70))
    self.nevents += 1
    leptonic = False
    tauonic = False
    bosons = [ ]
    taus = [ ]
    particles = Collection(event,'GenPart')
    #particles = Collection(event,'LHEPart')
    print(" \033[4m%7s %7s %7s %7s %7s %7s %7s %7s %8s %9s %10s  \033[0m"%(
      "index","pdgId","moth","mothId","dR","pt","eta","status","prompt","taudecay","last copy"))
    for i, particle in enumerate(particles):
      mothidx  = particle.genPartIdxMother
      eta      = max(-999,min(999,particle.eta))
      prompt   = particle.statusflag('isPrompt') #hasbit(particle.statusFlags,0)
      taudecay = particle.statusflag('isTauDecayProduct') #hasbit(particle.statusFlags,2)
      lastcopy = particle.statusflag('isLastCopy') #hasbit(particle.statusFlags,13)
      if 0<=mothidx<len(particles):
        moth    = particles[mothidx]
        mothpid = moth.pdgId
        mothdR  = max(-999,min(999,particle.DeltaR(moth))) #particle.p4().DeltaR(moth.p4())
        print(" %7d %7d %7d %7d %7.2f %7.2f %7.2f %7d %8s %9s %10s"%(
          i,particle.pdgId,mothidx,mothpid,mothdR,particle.pt,eta,particle.status,prompt,taudecay,lastcopy))
      else:
        print(" %7d %7d %7s %7s %7s %7.2f %7.2f %7d %8s %9s %10s"%(
          i,particle.pdgId,"","","",particle.pt,eta,particle.status,prompt,taudecay,lastcopy))
      if lastcopy:
        if abs(particle.pdgId) in [11,13,15]:
          leptonic = True
          if abs(particle.pdgId)==15:
            tauonic = True
            taus.append(particle)
        elif abs(particle.pdgId) in [23,24]:
          bosons.append(particle)
    for boson in bosons: # print production chain
      print("Boson production:")
      chain = getprodchain(boson,particles)[:-3]
      print(chain+getdecaychain(boson,particles,indent=len(chain),depth=0))
    for tau in taus: # print decay chain
      print("Tau decay:")
      print(getdecaychain(tau,particles))
    if leptonic:
      self.nleptonic += 1
    if tauonic:
      self.ntauonic += 1
  
  def endJob(self):
    print('\n'+'–'*96)
    if self.nevents>0:
      print("  %-10s %4d / %-4d (%.1f%%)"%('Tauonic: ',self.ntauonic, self.nevents,100.0*self.ntauonic/self.nevents))
      print("  %-10s %4d / %-4d (%.1f%%)"%('Leptonic:',self.nleptonic,self.nevents,100.0*self.nleptonic/self.nevents))
    print("%s Done %s\n"%('–'*10,'–'*80))
  

# PROCESS NANOAOD
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('-i', '--infiles', nargs='+')
parser.add_argument('-o', '--outdir', default='.')
parser.add_argument('-n', '--maxevts', type=int, default=20)
args = parser.parse_args()
infiles = args.infiles or [
  'root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18NanoAODv9/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/280000/525CD279-3344-6043-98B9-2EA8A96623E4.root',
]
processor = PostProcessor(args.outdir,infiles,noOut=True,modules=[LHEDumper()],maxEntries=args.maxevts)
processor.run()
