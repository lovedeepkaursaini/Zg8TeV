import FWCore.ParameterSet.Config as cms

from RecoLocalCalo.Configuration.ecalLocalRecoSequence_cff import *

process = cms.Process("ggKIT")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100)
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.GlobalTag.globaltag = cms.string('FT_53_V21_AN6::All')

#jec from sqlite
process.load("CondCore.DBCommon.CondDBCommon_cfi")
from CondCore.DBCommon.CondDBSetup_cfi import *
process.jec = cms.ESSource("PoolDBESSource",CondDBSetup,
			connect = cms.string('sqlite:Winter14_V8_DATA.db'),
			toGet = cms.VPSet(
				cms.PSet(
					record = cms.string('JetCorrectionsRecord'),
					tag = cms.string('JetCorrectorParametersCollection_Winter14_V8_DATA_AK5PFchs'),
					label = cms.untracked.string('AK4PFchs')
					),
				cms.PSet(
					record = cms.string('JetCorrectionsRecord'),
					tag = cms.string('JetCorrectorParametersCollection_Winter14_V8_DATA_AK7PFchs'),
					label = cms.untracked.string('AK8PFchs')
					)
				)
			)
process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    )
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
    #'/store/data/Run2012D/DoubleElectron/AOD/22Jan2013-v1/10027/FA331ACE-0B90-E211-9FE3-00261894393E.root'
    'file:/uscmst1b_scratch/lpc1/old_scratch/lpceg/yurii/lovedeep/ggNtuples/Nov26ggNtuples/CMSSW_5_3_12/src/ggAnalysis/ggNtuplizer/test/007A5A14-1069-E211-8BDD-0025905964B4.root'
    ), 
                            noEventSort = cms.untracked.bool(True),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
                            )

process.load("PhysicsTools.PatAlgos.patSequences_cff")
# Configure PAT to use PF2PAT instead of AOD sources
# this function will modify the PAT sequences.


# Trigger matching
process.load("ggAnalysis.ggNtuplizer.ggPatTriggerMatching_cff")

from PhysicsTools.PatAlgos.tools.cmsswVersionTools import *
#runOnData(process, ['All'], outputInProcess = False)

# Jets
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
process.patJetCorrFactors.useRho = cms.bool(True)

process.ak5PFJets.doAreaFastjet = True
process.patJets.addTagInfos = cms.bool(True)

# Taus
from PhysicsTools.PatAlgos.tools.tauTools import *
process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
process.cleanPatTaus.preselection = cms.string(' tauID("decayModeFinding") > 0.5 ')
process.cleanPatTaus.finalCut     = cms.string(' pt > 15.0 & abs(eta) < 2.5 ')
process.load("ggAnalysis.ggNtuplizer.ggTau_cff")

#process.patJetCorrFactors.levels = ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']
#process.patJetCorrFactors.rho = cms.InputTag('kt6PFJets25','rho')
#-------------------------------------
## Output Module Configuration (expects a path 'p')
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
process.out = cms.OutputModule("PoolOutputModule",
  fileName = cms.untracked.string("output.root"),
# save only events passing the full path
  SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
# save PAT Layer 1 output; you need a '*' to
# unpack the list of commands 'patEventContent'
  outputCommands = cms.untracked.vstring('drop *', *patEventContent)
)
from PhysicsTools.PatAlgos.tools.pfTools import *

postfix = "PFlow"
jetAlgo="AK5"
usePF2PAT(process,runPF2PAT=True, jetAlgo=jetAlgo, runOnMC=False, postfix=postfix)
# Add PF2PAT output to the created file
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContentNoCleaning
#process.load("CommonTools.ParticleFlow.PF2PAT_EventContent_cff")
#process.out.outputCommands =  cms.untracked.vstring('drop *')
process.out.outputCommands = cms.untracked.vstring('drop *',
                                                   'keep recoPFCandidates_particleFlow_*_*',
                                                   *patEventContentNoCleaning )


# top projections in PF2PAT:
getattr(process,"pfNoPileUp"+postfix).enable = True
getattr(process,"pfNoMuon"+postfix).enable = True
getattr(process,"pfNoElectron"+postfix).enable = True
getattr(process,"pfNoTau"+postfix).enable = False
getattr(process,"pfNoJet"+postfix).enable = True

# verbose flags for the PF2PAT modules
getattr(process,"pfNoMuon"+postfix).verbose = False

# enable delta beta correction for muon selection in PF2PAT?
getattr(process,"pfIsolatedMuons"+postfix).doDeltaBetaCorrection = False
## Jet energy corrections
jetCorrectionsAK5 = ('AK5PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual'])
jetCorrectionsAK7 = ('AK7PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual'])

## b-tag infos
bTagInfos = ['impactParameterTagInfos','secondaryVertexTagInfos']
## b-tag discriminators
bTagDiscriminators = ['jetProbabilityBJetTags','combinedSecondaryVertexBJetTags']

from PhysicsTools.PatAlgos.tools.jetTools import *
addJetCollection(process,
                 cms.InputTag('ak5PFJets'),
                 'AK5', 'PF',
                 doJTA        = True,
                 doBTagging   = True,
                 jetCorrLabel = jetCorrectionsAK5,#('AK5PF', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual'])),
                 doType1MET   = False,
                 doL1Cleaning = False,
                 doL1Counters = False,
                 genJetCollection = cms.InputTag("ak5GenJets"),
                 doJetID      = False
                 )

# CA8 jets (Gen and Reco)
from RecoJets.JetProducers.ca4GenJets_cfi import ca4GenJets
process.ca8GenJetsNoNu = ca4GenJets.clone(
    rParam = cms.double(0.8),
    src = cms.InputTag("genParticlesForJetsNoNu")
)

from RecoJets.JetProducers.ca4PFJets_cfi import ca4PFJets
process.ca8PFJets = ca4PFJets.clone(
    rParam = cms.double(0.8),
    src = getattr(process,"pfJets"+postfix).src,
    srcPVs = getattr(process,"pfJets"+postfix).srcPVs,
    doAreaFastjet = cms.bool(True),
    jetPtMin = cms.double(150)
)

## CA8 pruned jets (Gen and Reco)
from RecoJets.JetProducers.ca4GenJets_cfi import ca4GenJets
from RecoJets.JetProducers.SubJetParameters_cfi import SubJetParameters
process.ca8GenJetsNoNuPruned = ca4GenJets.clone(
    SubJetParameters,
    rParam = cms.double(0.8),
    src = cms.InputTag("genParticlesForJetsNoNu"),
    usePruning = cms.bool(True),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets")
)
from RecoJets.JetProducers.ak5PFJetsPruned_cfi import ak5PFJetsPruned
process.ca8PFJetsPruned = ak5PFJetsPruned.clone(
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam = cms.double(0.8),
    src = getattr(process,"pfJets"+postfix).src,
    srcPVs = getattr(process,"pfJets"+postfix).srcPVs,
    doAreaFastjet = cms.bool(True),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets"),
    jetPtMin = cms.double(150.)
)

## PATify CA8 jets
switchJetCollection(process,
    cms.InputTag('ca8PFJets'),
    doJTA = True,
    doBTagging = True,
    btagInfo = bTagInfos,
    btagdiscriminators = bTagDiscriminators,
    jetCorrLabel = jetCorrectionsAK7,
    doType1MET = False,
    genJetCollection = cms.InputTag('ca8GenJetsNoNu'),
    doJetID = False
)
addJetCollection(
    process,
    cms.InputTag('ca8PFJets'),
    'CA8','PF',
    doJTA=False,
    doBTagging=True,
    btagInfo=bTagInfos,
    btagdiscriminators=bTagDiscriminators,
    jetCorrLabel=jetCorrectionsAK7,
    doType1MET=False,
    doL1Cleaning=False,
    doL1Counters=False,
    doJetID=False,
    genJetCollection=cms.InputTag('ca8GenJetsNoNu')
)
addJetCollection(
    process,
    cms.InputTag('ca8PFJetsPruned'),
    'CA8Pruned','PF',
    doJTA=False,
    doBTagging=True,
    btagInfo=bTagInfos,
    btagdiscriminators=bTagDiscriminators,
    jetCorrLabel=jetCorrectionsAK7,
    doType1MET=False,
    doL1Cleaning=False,
    doL1Counters=False,
    doJetID=False,
    genJetCollection=cms.InputTag('ca8GenJetsNoNu')
)
addJetCollection(
    process,
    cms.InputTag('ca8PFJetsPruned','SubJets'),
    'CA8PrunedSubJets', 'PF',
    doJTA=True,
    doBTagging=True,
    btagInfo=bTagInfos,
    btagdiscriminators=bTagDiscriminators,
    jetCorrLabel=jetCorrectionsAK5,
    doType1MET=False,
    doL1Cleaning=False,
    doL1Counters=False,
    doJetID=False,
    genJetCollection=cms.InputTag('ca8GenJetsNoNuPruned','SubJets')
)

process.jetSeq = cms.Sequence(
  process.ca8PFJets
  + process.ca8PFJetsPruned
)

## N-subjettiness
process.selectedPatJetsCA8withNsub = cms.EDProducer("NjettinessAdder",
                              src=cms.InputTag("selectedPatJetsCA8PF"),
                              cone=cms.double(0.8)
                              )

## Establish references between PAT fat jets and PAT subjets using the BoostedJetMerger
process.selectedPatJetsCA8PrunedPFPacked = cms.EDProducer("BoostedJetMerger",
  jetSrc=cms.InputTag("selectedPatJetsCA8PrunedPF"),
  subjetSrc=cms.InputTag("selectedPatJetsCA8PrunedSubJetsPF")
)

from PhysicsTools.PatAlgos.tools.coreTools import *
removeMCMatching(process, names=['All'], outputModules=[])

process.load("JetMETCorrections.Type1MET.pfMETCorrections_cff")
process.pfJetMETcorr.jetCorrLabel = cms.string("ak5PFL1FastL2L3Residual")

# load the coreTools of PAT
from PhysicsTools.PatAlgos.tools.metTools import *
addTcMET(process, 'TC')
addPfMET(process, 'PF')

process.load("JetMETCorrections.Type1MET.pfMETCorrectionType0_cfi")
process.pfType0pfcp1CorrectedMet = process.pfType1CorrectedMet.clone()
process.pfType0pfcp1CorrectedMet.srcType1Corrections = cms.VInputTag(
    cms.InputTag('pfMETcorrType0'),
    cms.InputTag('pfJetMETcorr', 'type1')
    )

process.patMETsType0pfcp1PF = process.patMETsPF.clone()
process.patMETsType0pfcp1PF.metSource = cms.InputTag("pfType0pfcp1CorrectedMet")

process.produceType0MET = cms.Sequence(
    process.pfType0pfcp1CorrectedMet*
    process.patMETsType0pfcp1PF
    )

#process.patJetGenJetMatch.matched = cms.InputTag('iterativeCone5GenJets')

process.cleanPatPhotons.checkOverlaps.electrons.requireNoOverlaps = cms.bool(False)

# PF isolations for electrons and muons
from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFElectronIso, setupPFPhotonIso
process.eleIsoSequence = setupPFElectronIso(process, 'gsfElectrons')
process.phoIsoSequence = setupPFPhotonIso(process, 'photons')

process.load("ggAnalysis.ggNtuplizer.ggNtuplizer_crab_cfi")
process.ggNtuplizer.doGenParticles = cms.bool(False)
process.ggNtuplizer.jetSrc = cms.InputTag("selectedPatJetsAK5PF")
process.ggNtuplizer.triggerResults = cms.InputTag("TriggerResults::HLT")
process.ggNtuplizer.getBlocks=cms.bool(False)
process.ggNtuplizer.useAllPF=cms.bool(False)
process.ggNtuplizer.dumpTrks=cms.bool(False)
process.ggNtuplizer.dumpSubJets=cms.bool(True)
process.TFileService = cms.Service("TFileService", fileName = cms.string('ggtree_data.root'))

# electron energy regression
process.load("EgammaAnalysis.ElectronTools.electronRegressionEnergyProducer_cfi")
process.eleRegressionEnergy.energyRegressionType = cms.uint32(2)

process.load("Configuration.StandardSequences.Services_cff")
process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
                                                   calibratedPatElectrons = cms.PSet(
    initialSeed = cms.untracked.uint32(1),
    engineName = cms.untracked.string('TRandom3')
    ),
                                                   )

process.patElectrons.electronIDSources = cms.PSet(
    mvaTrigV0 = cms.InputTag("mvaTrigV0"),
    mvaNonTrigV0 = cms.InputTag("mvaNonTrigV0")
    )

process.load("EgammaAnalysis.ElectronTools.calibratedPatElectrons_cfi")
process.calibratedPatElectrons.isMC = cms.bool(False)
process.calibratedPatElectrons.inputDataset = cms.string("22Jan2013ReReco")
process.calibratedPatElectrons.correctionsType = cms.int32(2)
process.calibratedPatElectrons.combinationType = cms.int32(3)

process.load("ggAnalysis.ggNtuplizer.ggRhoFastJet_cff")
#process.load("ggAnalysis.ggNtuplizer.ggMergedJets_data_cff")
process.load("ggAnalysis.ggNtuplizer.ggEleID_cff")
process.load("ggAnalysis.ggNtuplizer.ggMETFilters_cff")
process.load("ggAnalysis.ggNtuplizer.ggBoostedEleModIso_cff")

process.patElectrons.userIsolation.user = cms.VPSet(
    cms.PSet(src = cms.InputTag("modElectronIso","track")),
    cms.PSet(src = cms.InputTag("modElectronIso","ecal")),
    cms.PSet(src = cms.InputTag("modElectronIso","hcalDepth1"))
    )

process.patElectrons.electronIDSources = cms.PSet(
    mvaTrigV0 = cms.InputTag("mvaTrigV0"),
    mvaNonTrigV0 = cms.InputTag("mvaNonTrigV0")
    )


# Trigger requirements
process.load('HLTrigger.HLTfilters.hltHighLevel_cfi')
process.leptonHLTFilter = copy.deepcopy(process.hltHighLevel)
process.leptonHLTFilter.throw = cms.bool(False)
process.leptonHLTFilter.HLTPaths = ['HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v*','HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v*','HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*','HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*','HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass50_v*','HLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50_v*','HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_v*','HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_v*','HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v*','HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_v*','HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v*']

process.scrapingVeto = cms.EDFilter("FilterOutScraping",
                                             applyfilter = cms.untracked.bool(True),
                                             debugOn = cms.untracked.bool(False),
                                             numtrack = cms.untracked.uint32(10),
                                             thresh = cms.untracked.double(0.25)
                                             )
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                      vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                      minimumNDOF = cms.uint32(4) ,
                      maxAbsZ = cms.double(24),
                      maxd0 = cms.double(2))
process.load('CommonTools/RecoAlgos/HBHENoiseFilter_cfi')
process.load("RecoMET.METFilters.hcalLaserEventFilter_cfi")
process.hcalLaserEventFilter.vetoByRunEventNumber=cms.untracked.bool(False)
process.hcalLaserEventFilter.vetoByHBHEOccupancy=cms.untracked.bool(True)
process.load('RecoMET.METFilters.eeBadScFilter_cfi')
process.load('RecoMET.METAnalyzers.CSCHaloFilter_cfi')
process.load('RecoMET.METFilters.EcalDeadCellTriggerPrimitiveFilter_cfi')
process.EcalDeadCellTriggerPrimitiveFilter.tpDigiCollection = cms.InputTag("ecalTPSkimNA")
process.filterSequence = cms.Sequence(
              process.scrapingVeto *
              process.primaryVertexFilter*
              process.HBHENoiseFilter*
              process.hcalLaserEventFilter*
              process.CSCTightHaloFilter*
              process.eeBadScFilter*
              process.EcalDeadCellTriggerPrimitiveFilter
          )


#------------------------------------- 
process.p = cms.Path(
    #process.leptonHLTFilter*
    process.filterSequence*
    process.fjSequence*
    process.ak5PFJets*
    process.pfNoPileUpSequence* ###########
    process.pfParticleSelectionSequence*
    process.ggBoostedEleModIsoSequence*
    process.eleMVAID*
    process.type0PFMEtCorrection*
    getattr(process,"patPF2PATSequence"+postfix)*
    process.jetSeq*
    getattr(process,"patDefaultSequence")*
    process.selectedPatJetsCA8PrunedPFPacked*
    process.selectedPatJetsCA8withNsub*
   #  process.patDefaultSequence*
    #process.Njettiness*
    process.produceType0MET*
    process.eleIsoSequence*
    process.phoIsoSequence*
#    process.ca8Jets* ###########
#    process.QuarkGluonTagger*
    process.eleRegressionEnergy*
    process.calibratedPatElectrons*
    process.ggTriggerSequence*
    process.ggMETFiltersSequence*
    process.recoTauClassicHPSSequence*
    process.ggNtuplizer)

del process.out
#process.out_step = cms.EndPath(process.output)
