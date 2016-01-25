import FWCore.ParameterSet.Config as cms

process = cms.Process('GEN')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('Hadronizer_Tune4C_8TeV_aMCatNLO_LHE_pythia8_cff_Z')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic7TeVCollision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.EventContent.EventContent_cff')

process.maxEvents = cms.untracked.PSet(
  input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("LHESource",
                            fileNames = cms.untracked.vstring('file:cmsgrid_final.lhe')
                            )

process.options = cms.untracked.PSet(
)
process.genskim = cms.EDFilter("CandViewSelector",
                               src = cms.InputTag("genParticles"),
                               cut = cms.string("pt() > 100 & (abs(pdgId) = 22 || abs(pdgId) = 23 || abs(pdgId) = 24)")
                               )
process.genfilter = cms.EDFilter("CandViewCountFilter",
                                 src = cms.InputTag("genskim"),
                                 minNumber = cms.uint32(1)
                                 )

# Production Info
process.configurationMetadata = cms.untracked.PSet(
  version = cms.untracked.string('$Revision: 1.0 $'),
  annotation = cms.untracked.string('Hadronizer_TuneCUETP8M1_13TeV_MLM_4f_max1j_LHE_pythia8_cff. nevts:-1'),
  name = cms.untracked.string('Tuning')
)

# Output definition

process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
  splitLevel = cms.untracked.int32(0),
  eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
  outputCommands = process.RAWSIMEventContent.outputCommands,
  fileName = cms.untracked.string('test.root'),
  dataset = cms.untracked.PSet(
    filterName = cms.untracked.string(''),
    dataTier = cms.untracked.string('GEN')
  ),
  SelectEvents = cms.untracked.PSet(
    SelectEvents = cms.vstring('generation_step')
  )
)

# Additional output definition

# Other statements
#process.GlobalTag.globaltag = 'START311_V2::All'

# Path and EndPath definitions
process.generation_step = cms.Path(process.generator*process.pgen*process.genskim*process.genfilter)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step = cms.EndPath(process.RAWSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.endjob_step,process.RAWSIMoutput_step)
