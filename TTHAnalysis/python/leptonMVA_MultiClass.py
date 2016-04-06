from array import array
from math import *
from CMGTools.TTHAnalysis.analyzers.ntupleTypes import ptRelv2, jetLepAwareJEC
from PhysicsTools.HeppyCore.utils.deltar import deltaR
from CMGTools.TTHAnalysis.signedSip import qualityTrk

import ROOT
#import os
#if "/smearer_cc.so" not in ROOT.gSystem.GetLibraries(): 
#    ROOT.gROOT.ProcessLine(".L %s/src/CMGTools/TTHAnalysis/python/plotter/smearer.cc+" % os.environ['CMSSW_BASE']);
#if "/mcCorrections_cc.so" not in ROOT.gSystem.GetLibraries(): 
#    ROOT.gROOT.ProcessLine(".L %s/src/CMGTools/TTHAnalysis/python/plotter/mcCorrections.cc+" % os.environ['CMSSW_BASE']);


class MVAVar:
    def __init__(self,name,func,corrfunc=None):
        self.name = name
        self.var  = array('f',[0.])
        self.func = func
        self.corrfunc = corrfunc
    def set(self,lep,ncorr): ## apply correction ncorr times
        self.var[0] = self.func(lep)
        if self.corrfunc:
            for i in range(ncorr):
                self.var[0] = self.corrfunc(self.var[0],lep) # pass lambda function as corrfunc if needed
class MVATool:
    def __init__(self,name,xml,specs,vars,nclasses=1):
        self.name = name
        self.reader = ROOT.TMVA.Reader("Silent")
        self.specs = specs
        self.vars  = vars
        self.nClasses = nClasses
        for s in specs: self.reader.AddSpectator(s.name,s.var)
        for v in vars:  self.reader.AddVariable(v.name,v.var)
        #print "Would like to load %s from %s! " % (name,xml)
        self.reader.BookMVA(name,xml)
    def __call__(self,lep,ncorr): ## apply correction ncorr times
        for s in self.specs: s.set(lep,ncorr)
        for s in self.vars:  s.set(lep,ncorr)
        return self.reader.EvaluateMVA(self.name) if self.nClasses==1 else self.reader.EvaluateMulticlass(self.name) # returns vector if multiclass
class CategorizedMVA:
    def __init__(self,catMvaPairs):
        self.catMvaPairs = catMvaPairs
        if len(self.catMvaPairs)<1: raise RuntimeError, 'No categories passed to CategorizedMVA'
        self.nClasses = self.catMvaPairs[0][1].nClasses
        if any([(x[1].nClasses!=self.nClasses) for x in self.catMvaPairs]): raise RuntimeError, 'Wrong multiclass configuration in CategorizedMVA'
    def __call__(self,lep,ncorr):
        for c,m in self.catMvaPairs:
            if c(lep): return m(lep,ncorr)
        return -99. if self.nClasses==1 else [-99.]*self.nClasses

_CommonSpect = [ 
]
_CommonVars = {
'forMoriond16':[ 
    MVAVar("LepGood_pt",lambda x: x.pt()),
    MVAVar("LepGood_eta",lambda x: x.eta()),
    MVAVar("LepGood_jetNDauChargedMVASel",lambda lepton: sum((deltaR(x.eta(),x.phi(),lepton.jet.eta(),lepton.jet.phi())<=0.4 and x.charge()!=0 and x.fromPV()>1 and qualityTrk(x.pseudoTrack(),lepton.associatedVertex)) for x in lepton.jet.daughterPtrVector()) if hasattr(lepton,'jet') and lepton.jet != lepton else 0),
    MVAVar("LepGood_miniRelIsoCharged",lambda x: getattr(x,'miniAbsIsoCharged',-99)/x.pt()), 
    MVAVar("LepGood_miniRelIsoNeutral",lambda x: getattr(x,'miniAbsIsoNeutral',-99)/x.pt()), 
    MVAVar("LepGood_jetPtRelv2", lambda x : ptRelv2(x) if hasattr(x,'jet') else -1),
    MVAVar("LepGood_jetPtRatio := min(LepGood_jetPtRatiov2,1.5)", lambda x : min((x.pt()/jetLepAwareJEC(x).Pt() if hasattr(x,'jet') else -1), 1.5)),
    MVAVar("LepGood_jetBTagCSV := max(LepGood_jetBTagCSV,0)", lambda x : max( (x.jet.btag('pfCombinedInclusiveSecondaryVertexV2BJetTags') if hasattr(x.jet, 'btag') else -99) ,0.)),
    MVAVar("LepGood_sip3d",lambda x: x.sip3D()),
    MVAVar("LepGood_dxy := log(abs(LepGood_dxy))",lambda x: log(abs(x.dxy()))),
    MVAVar("LepGood_dz  := log(abs(LepGood_dz))", lambda x: log(abs(x.dz()))),
 ],
}

_MuonVars = {
'forMoriond16': [
    MVAVar("LepGood_segmentCompatibility",lambda x: x.segmentCompatibility()), 
 ],
}

_ElectronVars = {
'forMoriond16': [
    MVAVar("LepGood_mvaIdSpring15",lambda x: x.mvaRun2("NonTrigSpring15")),
 ]

}


class LeptonMVA:
    def __init__(self, kind, basepath, isMC):
        global _CommonVars, _CommonSpect, _ElectronVars
        #print "Creating LeptonMVA of kind %s, base path %s" % (kind, basepath)
        self._isMC = isMC
        self._kind = kind
        muVars = _CommonVars[kind] + _MuonVars[kind]
        elVars = _CommonVars[kind] + _ElectronVars[kind]
        if 'forMoriond16' in self._kind:
            self.mu = CategorizedMVA([
                    ( lambda x: True, MVATool("BDTG",basepath%"mu",_CommonSpect,muVars) ),
                    ])
            self.el = CategorizedMVA([
                    ( lambda x: True, MVATool("BDTG",basepath%"el",_CommonSpect,elVars) ),
                    ])            
    def __call__(self,lep,ncorr="auto"):
        if ncorr == "auto": ncorr = 0 # (1 if self._isMC else 0)
        if   abs(lep.pdgId()) == 11: return self.el(lep,ncorr)
        elif abs(lep.pdgId()) == 13: return self.mu(lep,ncorr)
        else: return -99

