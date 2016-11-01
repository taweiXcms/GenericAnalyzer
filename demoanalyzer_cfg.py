import FWCore.ParameterSet.Config as cms
process = cms.Process("Demo")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#lxplus
'file:Pythia8_prompt_D0pt20p0_Pthat10_TuneCUETP8M1_8160GeV_cfi_evtgen130_py_GEN_SIM_PU.root'
    )
# skipEvents = cms.untracked.uint32(39)
)

process.genParticlePlusGEANT = cms.EDProducer("GenPlusSimParticleProducer",
        src           = cms.InputTag("g4SimHits"), # use "famosSimHits" for FAMOS
        setStatus     = cms.int32(8),             # set status = 8 for GEANT GPs
        filter        = cms.vstring("pt > 0.0"),  # just for testing (optional)
	    genParticles   = cms.InputTag("genParticles") # original genParticle list
	    #genParticles   = cms.InputTag("hiGenParticles") # original genParticle list
)

#process.demo = cms.EDAnalyzer('DemoAnalyzer')
process.demo = cms.EDAnalyzer('GenericAnalyzer',

    GenLabel        = cms.InputTag('genParticles'),
    #GenLabel        = cms.InputTag('hiGenParticles'),
    #GenLabel        = cms.InputTag('genParticlePlusGEANT'),
    MuonLabel       = cms.InputTag('muons'),
    #TrackLabel      = cms.InputTag('generalTracks'),
    TrackLabel      = cms.InputTag('hiGeneralTracks'),
    JetLabel      = cms.InputTag('ak3HiGenJets'),

	doHepMC = cms.untracked.bool(False),
	doGenParticle = cms.untracked.bool(True),
	doRecoMuon = cms.untracked.bool(False),
	#doRecoMuon = cms.untracked.bool(True)
	doRecoTrk = cms.untracked.bool(False),
	#doRecoTrk = cms.untracked.bool(True),
	doGenJet = cms.untracked.bool(False),

    PIDtoPrint = cms.vint32(421),
)
process.TFileService = cms.Service("TFileService",
      fileName = cms.string('results.root')      )

process.gp = cms.Path(process.genParticlePlusGEANT)
process.p = cms.Path(process.demo)
#process.schedule = cms.Schedule(process.gp, process.p)
process.schedule = cms.Schedule(process.p)
