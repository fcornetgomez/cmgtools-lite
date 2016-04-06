#!/usr/bin/env python
from CMGTools.TTHAnalysis.treeReAnalyzer import *
from array import array
from glob import glob
import os.path

import ROOT
from CMGTools.TTHAnalysis.leptonMVA import MVAVar, MVATool, CategorizedMVA

_CommonSpect = [ 
]
_CommonVars = {
'forMoriond16': [
    MVAVar("LepGood_pt",lambda x: x.pt),
    MVAVar("LepGood_eta",lambda x: x.eta),
    MVAVar("LepGood_jetNDauChargedMVASel",lambda x: x.jetNDauChargedMVASel),
    MVAVar("LepGood_miniRelIsoCharged",lambda x: x.miniRelIsoCharged),
    MVAVar("LepGood_miniRelIsoNeutral",lambda x: x.miniRelIsoNeutral),
    MVAVar("LepGood_jetPtRelv2",lambda x: x.jetPtRelv2),
    MVAVar("LepGood_jetPtRatio := min(LepGood_jetPtRatiov2,1.5)", lambda x : min(x.jetPtRatiov2,1.5)),
    MVAVar("LepGood_jetBTagCSV := max(LepGood_jetBTagCSV,0)", lambda x : max(x.jetBTagCSV,0.)),
    MVAVar("LepGood_sip3d",lambda x: x.sip3d),
    MVAVar("LepGood_dxy := log(abs(LepGood_dxy))",lambda x: log(abs(x.dxy))),
    MVAVar("LepGood_dz  := log(abs(LepGood_dz))", lambda x: log(abs(x.dz))),
]
}

_ElectronVars = {
 'forMoriond16': [
    MVAVar("LepGood_mvaIdSpring15",lambda x: x.mvaIdSpring15)
 ]
}

_MuonVars = [
    MVAVar("LepGood_segmentCompatibility",lambda x: x.segmentCompatibility)
]

class LeptonMVA:
    def __init__(self,basepath,training="forMoriond16",nClasses=1):
        global _CommonVars, _CommonSpect, _ElectronVars, _MuonVars, _SVVars
        if type(basepath) == tuple: basepathmu, basepathel  = basepath
        else:                       basepathmu, basepathel  = basepath, basepath
        print "Booking %s %s" % (training, basepath)
        if "forMoriond16" in training:
            muVars = _CommonVars[training][:] + _MuonVars[:]
            elVars = _CommonVars[training][:] + _ElectronVars[training][:]
        if not muVars:
            self.mu = lambda mu, ncorr : -37.0;
        elif "forMoriond16" in training:
            self.mu = CategorizedMVA([
                ( lambda x: True , MVATool("BDTG",basepathmu%"mu",_CommonSpect,muVars,nClasses) ),
            ])
        if not elVars:
            self.el = lambda el, ncorr : -37.0;
        elif "forMoriond16" in training:
            self.el = CategorizedMVA([
                ( lambda x: True, MVATool("BDTG",basepathel%"el",_CommonSpect,elVars,nClasses) ),
            ])
    def __call__(self,lep,ncorr=0):
        if   abs(lep.pdgId) == 11: return self.el(lep,ncorr)
        elif abs(lep.pdgId) == 13: return self.mu(lep,ncorr)
        else: return -99 if self.nClasses==1 else [-99]*self.nClasses

class LepMVAFriend:
    def __init__(self,path,training="forMoriond16",label="",nClasses=1):
        self.mva = LeptonMVA(path+"/%s_BDTG.weights.xml" if type(path) == str else path, training=training, nClasses=nClasses)
        self.label = label
    def listBranches(self):
        mylist = [ ("nLepGood","I") ]
        if self.nClasses>1:
            for i in xrange(self.nClasses): mylist.append(("LepGood_mva"+self.label+"_cl%d"%i,"F",8,"nLepGood"))
        else: mylist.append(("LepGood_mva"+self.label,"F",8,"nLepGood"))
        return mylist
    def __call__(self,event,ncorr="auto"):
        if ncorr == "auto": ncorr = 0
        lep = Collection(event,"LepGood","nLepGood",8)
        ret = { 'nLepGood' : event.nLepGood }
        if self.nClasses>1:
            allmvas = [self.mva(l, ncorr=0) for l in lep]
            for i in xrange(self.nClasses): ret['LepGood_mva'+self.label+'_cl%d'%i] = [ x[i] for x in allmvas ]
        else: ret['LepGood_mva'+self.label] = [ self.mva(l, ncorr=0) for l in lep ]
        return ret

if __name__ == '__main__':
    from sys import argv
    file = ROOT.TFile(argv[1])
    tree = file.Get("tree")
#    if len(argv) > 2: tree.AddFriend("sf/t",argv[2])
    tree.vectorTree = True
    class Tester(Module):
        def __init__(self, name, trees="new"):
            Module.__init__(self,name,None)
            self.mvas = {
                'forMoriond16' : LepMVAFriend(("/afs/cern.ch/user/p/peruzzi/work/cmgtools/CMSSW_7_4_14/src/CMGTools/TTHAnalysis/macros/leptons/weights/forMoriond16%s_BDTG.weights.xml",
                                               "/afs/cern.ch/user/p/peruzzi/work/cmgtools/CMSSW_7_4_14/src/CMGTools/TTHAnalysis/macros/leptons/weights/forMoriond16%s_BDTG.weights.xml",),
                                              training="forMoriond16",nClasses=1),
            }
        def analyze(self,ev):
            print "\nrun %6d lumi %4d event %d: leps %d" % (ev.run, ev.lumi, ev.evt, ev.nLepGood)
            lep = Collection(ev,"LepGood","nLepGood",8)
            for l,m in self.mvas.iteritems():
                print "%-10s: %s %s" % (l, m(ev), [ x.mvaTTH for x in lep ] )
    el = EventLoop([ Tester("tester", "new") ])
    el.loop([tree], maxEvents = 50)

        
