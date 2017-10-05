import FWCore.ParameterSet.Config as cms
process = cms.Process("Demo")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) )
process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#lxplus
#'file:/afs/cern.ch/work/t/twang/MITHIG/GenHIBmeson_20131220/BoostGentest/CMSSW_5_3_8_HI_patch2/src/test/localRun/RunMore/PyquenMix_embedHIJING_Bp2JpsiKp_5TeV.root'
#'file:/afs/cern.ch/work/t/twang/MITHIG/GenHIBmeson_20131220/BoostGentest_20140212/localRun/RunMore/PyquenMix_embedHIJING_Bp2JpsiKp_5TeV_boostedMC.root'
#'file:/afs/cern.ch/work/t/twang/MITHIG/GenHIBmeson_20131220/BoostGen5GeVB_20140214/localRun/RunMore/PyquenMix_embedHIJING_Bp2JpsiKp_Bpt5_5TeV_boostedMC.root'

#'file:/afs/cern.ch/work/t/twang/MITHIG/GenHIBmeson_20131220/BoostGen5GeVB_20140214/subGENSIM_20140219/subBu/step4/step4_RAW2DIGI_L1Reco_RECO.root'
#'file:/afs/cern.ch/work/t/twang/MITHIG/GenHIBmeson_20131220/BoostGen5GeVB_20140214/subGENSIM_20140219/subBdKs/HIJINGemb_BdJpsiKs_TuneZ2star_5TeV_cff_GEN_SIM.root'
#'file:/afs/cern.ch/work/t/twang/MITHIG/GenHIBmeson_20131220/BoostGen5GeVB_20140214/subGENSIM_20140219/subBdKs/HIJINGemb_BdJpsiKs_TuneZ2star_5TeV_cff_GEN_SIM_B0ToJpsiKs.root'
#'file:/afs/cern.ch/work/t/twang/MITHIG/GenHIBmeson_20131220/BoostGen5GeVB_20140214/subGENSIM_20140219/subBdKs/HIJINGemb_BdJpsiKs_TuneZ2star_5TeV_cff_GEN_SIM_Bd_JpsiKs_mumu.root'
#'file:/afs/cern.ch/work/t/twang/MITHIG/GenHIBmeson_20131220/BoostGen5GeVB_20140214/subGENSIM_20140219/subBdKs/step4/HIJINGemb_BdJpsiKs_TuneZ2star_5TeV_cff_step4_RAW2DIGI_L1Reco_RECO.root'

#'file:/afs/cern.ch/work/t/twang/MITHIG/GenHIBmeson_20131220/BoostGen5GeVB_20140214/PYTHIA_sub/B2Kstar/PYTHIA6_BdJpsiKstar_TuneZ2star_2760GeV_GENSIM.root'
#'file:/afs/cern.ch/work/t/twang/MITHIG/GenHIBmeson_20131220/BoostGen5GeVB_20140214/PYTHIA_sub/B2Phi/PYTHIA6_BsJpsiPhi_TuneZ2star_2760GeV_GENSIM.root'

#'file:/afs/cern.ch/user/t/twang/work/MITHIG/GenHIBmeson_20131220/PbPbB_CMSSW75X_20151026/ourCustomiizeFilter/CMSSW_7_5_3_patch1/src/test/Pyquen_BuKp_Unquenched_5020GeV_cfi_py_GEN_SIM_PU.root'
#'file:/afs/cern.ch/user/t/twang/work/MITHIG/GenHIBmeson_20131220/PbPbB_CMSSW75X_20151026/useOfficialFilter/CMSSW_7_5_3_patch1/src/test/Pyquen_BuKp_Unquenched_5020GeV_cfi_py_GEN_SIM_PU.root'
#'file:/afs/cern.ch/user/t/twang/work/MITHIG/GenHIBmeson_20131220/PbPbB_CMSSW75X_20151026/useOfficialFilter/CMSSW_7_5_3_patch1/src/test/Pyquen_BdKstar_Unquenched_5020GeV_cfi_py_GEN_SIM_PU.root'
#'file:/afs/cern.ch/user/t/twang/work/MITHIG/GenHIBmeson_20131220/PbPbB_CMSSW75X_20151026/useOfficialFilter/CMSSW_7_5_3_patch1/src/test/Pyquen_BsPhi_Unquenched_5020GeV_cfi_py_GEN_SIM_PU.root'

#'file:/afs/cern.ch/user/t/twang/work/MITHIG/GenHIBmeson_20131220/PbPbB_CMSSW75X_20151026/useOfficialFilter/CMSSW_7_5_3_patch1/src/test/Pythia8_BuToJpsiK_TuneCUEP8M1_5020GeV_BPHMod_filter_cfi_py_GEN_SIM_PU.root'
#'file:/afs/cern.ch/user/t/twang/work/MITHIG/GenSpace/pp5022017_20171001/CMSSW_9_2_12_patch1/src/test/Pythia8_prompt_D0pt0p0_Pthat100_TuneCUETP8M1_5020GeV_cfi_evtgen130_py_GEN_SIM.root'
'file:/afs/cern.ch/user/t/twang/work/MITHIG/GenSpace/pp5022017_20171001/CMSSW_9_2_12_patch1/src/test/Pythia8_prompt_Dspt0p0_Pthat0_TuneCUETP8M1_5020GeV_cfi_evtgen130_py_GEN_SIM.root'
#'file:/afs/cern.ch/user/t/twang/work/MITHIG/GenSpace/pp5022017_20171001/CMSSW_9_2_12_patch1/src/testphoton/Pythia8_prompt_D0pt0p0_Gammapt15p0_Pthat0_TuneCUETP8M1_5020GeV_cfi_evtgen130_py_GEN_SIM.root'
    )
# skipEvents = cms.untracked.uint32(39)
)

process.genParticlePlusGEANT = cms.EDProducer("GenPlusSimParticleProducer",
        src           = cms.InputTag("g4SimHits"), # use "famosSimHits" for FAMOS
        setStatus     = cms.int32(8),             # set status = 8 for GEANT GPs
        filter        = cms.vstring("pt > 0.0"),  # just for testing (optional)
	    genParticles   = cms.InputTag("genParticles") # original genParticle list
#	    genParticles   = cms.InputTag("hiGenParticles") # original genParticle list
)

#process.demo = cms.EDAnalyzer('DemoAnalyzer')
process.demo = cms.EDAnalyzer('GenericAnalyzer',

    GenLabel        = cms.InputTag('genParticles'),
#    GenLabel        = cms.InputTag('hiGenParticles'),
#    GenLabel        = cms.InputTag('genParticlePlusGEANT'),
    MuonLabel       = cms.InputTag('muons'),
#    TrackLabel      = cms.InputTag('generalTracks'),
    TrackLabel      = cms.InputTag('hiGeneralTracks'),
    JetLabel      = cms.InputTag('ak3HiGenJets'),

	doHepMC = cms.untracked.bool(False),
	doGenParticle = cms.untracked.bool(True),
	doGenJet = cms.untracked.bool(False),
	doRecoTrk = cms.untracked.bool(False),
#	doRecoTrk = cms.untracked.bool(True),
	doRecoMuon = cms.untracked.bool(False)
#	doRecoMuon = cms.untracked.bool(True)
)
process.TFileService = cms.Service("TFileService",
      fileName = cms.string('results.root')      )

process.gp = cms.Path(process.genParticlePlusGEANT)
process.p = cms.Path(process.demo)
#process.schedule = cms.Schedule(process.gp, process.p)
process.schedule = cms.Schedule(process.p)
