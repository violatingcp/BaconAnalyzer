import FWCore.ParameterSet.Config as cms

#from Configuration.Generator.Pythia8CommonSettings_cfi import *
#from Configuration.Generator.Pythia8CUEP8M1Settings_cfi import *
#from Configuration.Generator.Pythia8aMCatNLOSettings_cfi import *

pythia8CommonSettingsBlock = cms.PSet(
    pythia8CommonSettings = cms.vstring(
      'Main:timesAllowErrors = 10000',
      'Check:epTolErr = 0.01',
      'SLHA:keepSM = on',
      'SLHA:minMassSM = 1000.',
      'ParticleDecays:limitTau0 = on',
      'ParticleDecays:tau0Max = 10',
      'ParticleDecays:allowPhotonRadiation = on',
    )
)

pythia8CUEP8M1SettingsBlock = cms.PSet(
    pythia8CUEP8M1Settings = cms.vstring(
        'Tune:pp 14',
        'Tune:ee 7',
        'MultipartonInteractions:pT0Ref=2.4024',
        'MultipartonInteractions:ecmPow=0.25208',
        'MultipartonInteractions:expPow=1.6',
    )
)
pythia8aMCatNLOSettingsBlock = cms.PSet(
    pythia8aMCatNLOSettings = cms.vstring(
      'SpaceShower:pTmaxMatch = 1',
      'SpaceShower:pTmaxFudge = 1',
      'SpaceShower:MEcorrections = off',
      'TimeShower:pTmaxMatch = 1',
      'TimeShower:pTmaxFudge = 1',
      'TimeShower:MEcorrections = off',
      'TimeShower:globalRecoil = on',
      'TimeShower:limitPTmaxGlobal = on',
      'TimeShower:nMaxGlobalRecoil = 1',
      'TimeShower:globalRecoilMode = 2',
      'TimeShower:nMaxGlobalBranch = 1',
    )
)

generator = cms.EDFilter("Pythia8HadronizerFilter",
    maxEventsToPrint = cms.untracked.int32(1),
    pythiaPylistVerbosity = cms.untracked.int32(1),
    filterEfficiency = cms.untracked.double(1.0),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    comEnergy = cms.double(13000.),
    PythiaParameters = cms.PSet(
        pythia8CommonSettingsBlock,
        pythia8CUEP8M1SettingsBlock,
        pythia8aMCatNLOSettingsBlock,
        processParameters = cms.vstring(
            'JetMatching:setMad = off',
            'JetMatching:scheme = 1',
            'JetMatching:merge = on',
            'JetMatching:jetAlgorithm = 2',
            'JetMatching:etaJetMax = 999.',
            'JetMatching:coneRadius = 1.',
            'JetMatching:slowJetPower = 1',
            'JetMatching:qCut = 70.', #this is the actual merging scale
            'JetMatching:doFxFx = on',
            'JetMatching:qCutME = 10.',#this must match the ptj cut in the lhe generation step
            'JetMatching:nQmatch = 5', #4 corresponds to 4-flavour scheme (no matching of b-quarks), 5 for 5-flavour scheme
            'JetMatching:nPartonsNow = 0',
            'TimeShower:nPartonsInBorn = 0', #number of coloured particles (before resonance decays) in highest multiplicity born matrix element
            'JetMatching:nJetMax = 2', #number of partons in born matrix element for highest multiplicity
            ),
        parameterSets = cms.vstring('pythia8CommonSettings',
                                    'pythia8CUEP8M1Settings',
                                    'pythia8aMCatNLOSettings',
                                    'processParameters',
                                    )
        )
    )

