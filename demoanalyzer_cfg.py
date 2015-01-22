import FWCore.ParameterSet.Config as cms
process = cms.Process("Demo")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#lxplus
#'file:/afs/cern.ch/work/t/twang/MITHIG/GenHIBmeson_20131220/CMSSW_5_3_8_HI_patch2/src/test/crab_0_131222_235821/res/PyquenMix_embedHIJING_Bp2JpsiKp_5TeV_9_1_Sm0.root'
#'file:/afs/cern.ch/work/t/twang/MITHIG/GenHIBmeson_20131220/CMSSW_5_3_8_HI_patch2/src/test/localRun/PyquenMix_embedHIJING_Bp2JpsiKp_5TeV.root'
#'file:/afs/cern.ch/work/t/twang/MITHIG/GenHIBmeson_20131220/CMSSW_5_3_8_HI_patch2/src/test/localRun/PyquenMix_embedHIJING_Bd2JpsiKstar_5TeV.root'
#'file:/afs/cern.ch/work/t/twang/MITHIG/GenHIBmeson_20131220/CMSSW_5_3_8_HI_patch2/src/test/localRun/PyquenMix_embedHIJING_Bd2JpsiKs_5TeV.root'
#'file:/afs/cern.ch/work/t/twang/MITHIG/GenHIBmeson_20131220/CMSSW_5_3_8_HI_patch2/src/test/localRun/PyquenMix_embedHIJING_Bs2JpsiPhi_5TeV.root'
#'file:/afs/cern.ch/work/t/twang/MITHIG/GenHIBmeson_20131220/BoostGentest/CMSSW_5_3_8_HI_patch2/src/test/localRun/RunMore/PyquenMix_embedHIJING_Bp2JpsiKp_5TeV.root'
#'file:/afs/cern.ch/work/t/twang/MITHIG/GenHIBmeson_20131220/BoostGentest_20140212/localRun/RunMore/PyquenMix_embedHIJING_Bp2JpsiKp_5TeV_boostedMC.root'
#'file:/afs/cern.ch/work/t/twang/MITHIG/GenHIBmeson_20131220/BoostGen5GeVB_20140214/localRun/RunMore/PyquenMix_embedHIJING_Bp2JpsiKp_Bpt5_5TeV_boostedMC.root'
#'file:/afs/cern.ch/work/t/twang/MITHIG/GenHIBmeson_20131220/BoostGen5GeVB_20140214/BoostEmbCodes/subBu/HIJINGemb_BuJpsiK_TuneZ2star_5TeV_cff_GEN_SIM.root'
#'file:/afs/cern.ch/work/t/twang/MITHIG/GenHIBmeson_20131220/BoostGen5GeVB_20140214/BoostEmbCodes/subBd/HIJINGemb_BdJpsiKstar_TuneZ2star_5TeV_cff_GEN_SIM.root'
#'file:/afs/cern.ch/work/t/twang/MITHIG/GenHIBmeson_20131220/BoostGen5GeVB_20140214/BoostEmbCodes/subBs/HIJINGemb_BsJpsiPhi_TuneZ2star_5TeV_cff_GEN_SIM.root'
#'file:/afs/cern.ch/work/t/twang/MITHIG/GenHIBmeson_20131220/BoostGen5GeVB_20140214/subGENSIM_20140219/subBu/step4/step4_RAW2DIGI_L1Reco_RECO.root'
#'file:/afs/cern.ch/cms/Tutorials/TWIKI_DATA/RelValTTbar_RECO_413.root'
#'file:/afs/cern.ch/work/t/twang/MITHIG/GenHIBmeson_20131220/BoostGen5GeVB_20140214/subGENSIM_20140219/subBdKs/HIJINGemb_BdJpsiKs_TuneZ2star_5TeV_cff_GEN_SIM.root'
#'file:/afs/cern.ch/work/t/twang/MITHIG/GenHIBmeson_20131220/BoostGen5GeVB_20140214/subGENSIM_20140219/subBdKs/HIJINGemb_BdJpsiKs_TuneZ2star_5TeV_cff_GEN_SIM_B0ToJpsiKs.root'
#'file:/afs/cern.ch/work/t/twang/MITHIG/GenHIBmeson_20131220/BoostGen5GeVB_20140214/subGENSIM_20140219/subBdKs/HIJINGemb_BdJpsiKs_TuneZ2star_5TeV_cff_GEN_SIM_Bd_JpsiKs_mumu.root'
#'file:/afs/cern.ch/work/t/twang/MITHIG/GenHIBmeson_20131220/BoostGen5GeVB_20140214/subGENSIM_20140219/subBdKs/step4/HIJINGemb_BdJpsiKs_TuneZ2star_5TeV_cff_step4_RAW2DIGI_L1Reco_RECO.root'
#'file:/afs/cern.ch/work/t/twang/MITHIG/GenHIBmeson_20131220/BoostGen5GeVB_20140214/PYTHIA_sub/B2Kstar/PYTHIA6_BdJpsiKstar_TuneZ2star_2760GeV_GENSIM.root'
#'file:/afs/cern.ch/work/t/twang/MITHIG/GenHIBmeson_20131220/BoostGen5GeVB_20140214/PYTHIA_sub/B2Phi/PYTHIA6_BsJpsiPhi_TuneZ2star_2760GeV_GENSIM.root'
#cgate
#'file:/net/hisrv0001/home/tawei/twang/HItestSpace/0306_pp_Kstar/PYTHIA6_BdJpsiKstar_TuneZ2star_2760GeV_GENSIM_3_1_Itr.root'
#'file:/net/hisrv0001/home/tawei/twang/HItestSpace/0306_pp_Phi/PYTHIA6_BsJpsiPhi_TuneZ2star_2760GeV_GENSIM_1_1_0u6.root'
#'file:/net/hisrv0001/home/tawei/twang/HItestSpace/0306_pp_Phi/PYTHIA6_BsJpsiPhi_TuneZ2star_2760GeV_GENSIM_3_1_T6M.root'
'file:/net/hisrv0001/home/tawei/HeavyFlavor_20131030/Gen_HiBmeson_20150121/PbPbtest/CMSSW_5_3_20/src/test/RECO.root'
    )
# skipEvents = cms.untracked.uint32(39)
)

process.genParticlePlusGEANT = cms.EDProducer("GenPlusSimParticleProducer",
        src           = cms.InputTag("g4SimHits"), # use "famosSimHits" for FAMOS
        setStatus     = cms.int32(8),             # set status = 8 for GEANT GPs
        filter        = cms.vstring("pt > 0.0"),  # just for testing (optional)
#	    genParticles   = cms.InputTag("genParticles") # original genParticle list
	    genParticles   = cms.InputTag("hiGenParticles") # original genParticle list
)

#process.demo = cms.EDAnalyzer('DemoAnalyzer')
process.demo = cms.EDAnalyzer('GenericAnalyzer',

#    GenLabel        = cms.InputTag('genParticles'),
    GenLabel        = cms.InputTag('hiGenParticles'),
#    GenLabel        = cms.InputTag('genParticlePlusGEANT'),
    MuonLabel       = cms.InputTag('muons'),
#    TrackLabel      = cms.InputTag('generalTracks'),
    TrackLabel      = cms.InputTag('hiGeneralTracks'),
    JetLabel      = cms.InputTag('ak3HiGenJets'),

	doHepMC = cms.untracked.bool(False),
	doGenParticle = cms.untracked.bool(True),
	doGenJet = cms.untracked.bool(False),
#	doRecoTrk = cms.untracked.bool(False),
	doRecoTrk = cms.untracked.bool(True),
	doRecoMuon = cms.untracked.bool(False)
#	doRecoMuon = cms.untracked.bool(True)
)
process.TFileService = cms.Service("TFileService",
      fileName = cms.string('results.root')      )

process.gp = cms.Path(process.genParticlePlusGEANT)
process.p = cms.Path(process.demo)
#process.schedule = cms.Schedule(process.gp, process.p)
process.schedule = cms.Schedule(process.p)
