import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.categories.append('PATLayer0Summary')
process.MessageLogger.categories.append('PATSummaryTables')
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAny_cfi")
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAlong_cfi")
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorOpposite_cfi")
process.load("PhysicsTools.PatAlgos.patSequences_cff")

process.load("RecoMuon.DetLayers.muonDetLayerGeometry_cfi")


process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.GlobalTag.globaltag = cms.string('IDEAL_V9::All')
process.GlobalTag.globaltag = 'START52_V4::All'
process.load("Configuration.StandardSequences.MagneticField_cff")





process.MessageLogger = cms.Service("MessageLogger",
	destinations = cms.untracked.vstring("Log")
)


process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
#process.options = cms.untracked.PSet(SkipEvent = cms.untracked.vstring('ProductNotFound'))
fileList = cms.untracked.vstring()
#fileList.extend(['root://xrootd.unl.edu//store/relval/CMSSW_6_2_0_SLHC11/RelValSingleMuPt10/GEN-SIM-RECO/DES23_62_V1_UPG2023Muon-v1/00000/E8AD6977-6EC6-E311-9ED6-0025905A60EE.root'])
fileList.extend(['file:/net/scratch_cms/institut_3a/scheuch/SingleMuPt10to40_001.root'])
fileList.extend(['file:/net/scratch_cms/institut_3a/scheuch/SingleMuPt10to40_002.root'])
fileList.extend(['file:/net/scratch_cms/institut_3a/scheuch/SingleMuPt10to40_003.root'])
fileList.extend(['file:/net/scratch_cms/institut_3a/scheuch/SingleMuPt10to40_004.root'])
fileList.extend(['file:/net/scratch_cms/institut_3a/scheuch/SingleMuPt10to40_005.root'])
fileList.extend(['file:/net/scratch_cms/institut_3a/scheuch/SingleMuPt10to40_006.root'])
fileList.extend(['file:/net/scratch_cms/institut_3a/scheuch/SingleMuPt10to40_007.root'])
fileList.extend(['file:/net/scratch_cms/institut_3a/scheuch/SingleMuPt10to40_008.root'])
fileList.extend(['file:/net/scratch_cms/institut_3a/scheuch/SingleMuPt10to40_009.root'])
fileList.extend(['file:/net/scratch_cms/institut_3a/scheuch/SingleMuPt10to40_010.root'])
fileList.extend(['file:/net/scratch_cms/institut_3a/scheuch/SingleMuPt10to40_011.root'])
#fileList.extend(['file:Pt100_1.root'])
# fileList.extend(['file:Pt10_1.root'])
# fileList.extend(['file:Pt10_2.root'])
# fileList.extend(['file:Pt10_3.root'])
#fileList.extend(['file:Pt1_1.root'])
#fileList.extend(['file:Pt1_2.root'])
#fileList.extend(['file:Pt1_3.root'])

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = fileList
)

#process.load("L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskAlgoTrigConfig_cff")
#process.es_prefer_l1GtTriggerMaskAlgoTrig = cms.ESPrefer("L1GtTriggerMaskAlgoTrigTrivialProducer", "l1GtTriggerMaskAlgoTrig")

#process.load('HLTrigger.HLTfilters.hltLevel1GTSeed_cfi')
#process.hltLevel1GTSeed.L1TechTriggerSeeding = cms.bool(False)
#process.hltLevel1GTSeed.L1SeedsLogicalExpression = cms.string('L1_DoubleMu3')

process.load('PhysicsTools.PatAlgos.mcMatchLayer0.muonMatch_cfi')
process.muonMatch.src = cms.InputTag("l1extraParticles")
process.muonMatch.mcPdgId = cms.vint32(13)
process.muonMatch.checkCharge = cms.bool(False)
process.muonMatch.resolveAmbiguities = cms.bool(False)
process.muonMatch.maxDeltaR = cms.double(10.)
process.muonMatch.maxDPtRel = cms.double(100.)
process.muonMatch.resolveByMatchQuality = cms.bool(True)

#process.load("L1TriggerConfig.L1GtConfigProducers.l1GtTriggerMenuXml_cfi")
#process.l1GtTriggerMenuXml._errorstr = 0;

#process.load('L1Trigger.Skimmer.l1Filter_cfi')
#process.l1Filter.algorithms = cms.vstring('L1_DoubleMu3')

process.load("MuonAnalysis.MuonAssociators.muonL1Match_cfi")
#process.patDefaultSequence.replace(process.muonTrigMatchHLT1MuonNonIso, process.muonTrigMatchHLT1MuonNonIso + process.muonL1Match)
#process.allLayer1Muons.trigPrimMatch += [ cms.InputTag("muonL1Match") ]


#process.selectedMuonsGenParticlesMatchNew = cms.EDProducer("MCTruthDeltaRMatcherNew",
#                                              src = cms.InputTag("muons"),
#                                              matched = cms.InputTag("genParticles"),
#                                              distMin = cms.double(0.15),
#                                              matchPDGId = cms.vint32(13)
#)

process.demo = cms.EDAnalyzer('L1Match'
)

#process.demo2 = cms.EDAnalyzer('MuonMatchPre'
#)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('AnalysisResult.root')
)

process.p = cms.Path(process.muonMatch*process.demo) #process.l1GtTriggerMenuXml, process.l1Filter*  *process.hltLevel1GTSeed #process.muonL1Match*
