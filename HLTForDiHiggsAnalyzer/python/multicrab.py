from CRABClient.UserUtilities import config
config = config()

# General

config.General.transferLogs    = True
config.General.transferOutputs = True
config.General.workArea = 'multicrab_Efficiencies_ggHHTo4b_15Feb2023'

# JobType
config.JobType.pluginName  = 'Analysis'
config.JobType.psetName    = 'ConfFile_cfg.py' # Name of the CMSSW configuration file
config.JobType.numCores    = 1
config.JobType.disableAutomaticOutputCollection = True
config.JobType.outputFiles = ['efficiency.root']
config.JobType.inputFiles  = ['setup_dev_CMSSW_13_0_0_GRun_cff.py']
config.JobType.maxMemoryMB = 2500

# Data
#config.Data.inputDataset          = "/GluGlutoHHto4B_kl-0p00_kt-1p00_c2-0p00_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22EEMiniAODv3-Poisson70KeepRAW_124X_mcRun3_2022_realistic_postEE_v1-v2/MINIAODSIM"
#config.Data.secondaryInputDataset = "/GluGlutoHHto4B_kl-0p00_kt-1p00_c2-0p00_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22EEDR-Poisson70KeepRAW_124X_mcRun3_2022_realistic_postEE_v1-v2/GEN-SIM-RAW" 
#config.Data.outputDatasetTag      = "GluGlutoHHto4B_kl-0p00_kt-1p00_c2-0p00_TuneCP5_13p6TeV_powheg-pythia8"
#config.General.requestName        = "GluGlutoHHto4B_kl-0p00_kt-1p00_c2-0p00_TuneCP5_13p6TeV_powheg-pythia8"


#config.Data.secondaryInputDataset = "/GluGlutoHHto4B_kl-0p00_kt-1p00_c2-1p00_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22EEDR-Poisson70KeepRAW_124X_mcRun3_2022_realistic_postEE_v1-v2/GEN-SIM-RAW"
#config.Data.inputDataset          = "/GluGlutoHHto4B_kl-0p00_kt-1p00_c2-1p00_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22EEMiniAODv3-Poisson70KeepRAW_124X_mcRun3_2022_realistic_postEE_v1-v2/MINIAODSIM"
#config.Data.outputDatasetTag      = "GluGlutoHHto4B_kl-0p00_kt-1p00_c2-1p00_TuneCP5_13p6TeV_powheg-pythia8"
#config.General.requestName        = "GluGlutoHHto4B_kl-0p00_kt-1p00_c2-1p00_TuneCP5_13p6TeV_powheg-pythia8"


#config.Data.secondaryInputDataset = "/GluGlutoHHto4B_kl-1p00_kt-1p00_c2-0p00_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22EEDR-Poisson70KeepRAW_124X_mcRun3_2022_realistic_postEE_v1-v2/GEN-SIM-RAW"
#config.Data.inputDataset          = "/GluGlutoHHto4B_kl-1p00_kt-1p00_c2-0p00_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22EEMiniAODv3-Poisson70KeepRAW_124X_mcRun3_2022_realistic_postEE_v1-v2/MINIAODSIM"
#config.Data.outputDatasetTag      = "GluGlutoHHto4B_kl-1p00_kt-1p00_c2-0p00_TuneCP5_13p6TeV_powheg-pythia8"
#config.General.requestName        = "GluGlutoHHto4B_kl-1p00_kt-1p00_c2-0p00_TuneCP5_13p6TeV_powheg-pythia8"


#config.Data.secondaryInputDataset = "/GluGlutoHHto4B_kl-1p00_kt-1p00_c2-0p35_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22EEDR-Poisson70KeepRAW_124X_mcRun3_2022_realistic_postEE_v1-v2/GEN-SIM-RAW"
#config.Data.inputDataset          = "/GluGlutoHHto4B_kl-1p00_kt-1p00_c2-0p35_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22EEMiniAODv3-Poisson70KeepRAW_124X_mcRun3_2022_realistic_postEE_v1-v2/MINIAODSIM"
#config.Data.outputDatasetTag      = "GluGlutoHHto4B_kl-1p00_kt-1p00_c2-0p35_TuneCP5_13p6TeV_powheg-pythia8"
#config.General.requestName        = "GluGlutoHHto4B_kl-1p00_kt-1p00_c2-0p35_TuneCP5_13p6TeV_powheg-pythia8"


#config.Data.secondaryInputDataset = "/GluGlutoHHto4B_kl-1p00_kt-1p00_c2-3p00_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22EEDR-Poisson70KeepRAW_124X_mcRun3_2022_realistic_postEE_v1-v2/GEN-SIM-RAW"
#config.Data.inputDataset          = "/GluGlutoHHto4B_kl-1p00_kt-1p00_c2-3p00_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22EEMiniAODv3-Poisson70KeepRAW_124X_mcRun3_2022_realistic_postEE_v1-v2/MINIAODSIM"
#config.Data.outputDatasetTag      = "GluGlutoHHto4B_kl-1p00_kt-1p00_c2-3p00_TuneCP5_13p6TeV_powheg-pythia8"
#config.General.requestName        = "GluGlutoHHto4B_kl-1p00_kt-1p00_c2-3p00_TuneCP5_13p6TeV_powheg-pythia8"


#config.Data.secondaryInputDataset = "/GluGlutoHHto4B_kl-1p00_kt-1p00_c2-m2p00_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22EEDR-Poisson70KeepRAW_124X_mcRun3_2022_realistic_postEE_v1-v2/GEN-SIM-RAW"
#config.Data.inputDataset          = "/GluGlutoHHto4B_kl-1p00_kt-1p00_c2-m2p00_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22EEMiniAODv3-Poisson70KeepRAW_124X_mcRun3_2022_realistic_postEE_v1-v2/MINIAODSIM"
#config.Data.outputDatasetTag      = "GluGlutoHHto4B_kl-1p00_kt-1p00_c2-m2p00_TuneCP5_13p6TeV_powheg-pythia8"
#config.General.requestName        = "GluGlutoHHto4B_kl-1p00_kt-1p00_c2-m2p00_TuneCP5_13p6TeV_powheg-pythia8"


#config.Data.secondaryInputDataset = "/GluGlutoHHto4B_kl-2p45_kt-1p00_c2-0p00_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22EEDR-Poisson70KeepRAW_124X_mcRun3_2022_realistic_postEE_v1-v2/GEN-SIM-RAW"
#config.Data.inputDataset          = "/GluGlutoHHto4B_kl-2p45_kt-1p00_c2-0p00_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22EEMiniAODv3-Poisson70KeepRAW_124X_mcRun3_2022_realistic_postEE_v1-v2/MINIAODSIM"
#config.Data.outputDatasetTag      = "GluGlutoHHto4B_kl-2p45_kt-1p00_c2-0p00_TuneCP5_13p6TeV_powheg-pythia8"
#config.General.requestName        = "GluGlutoHHto4B_kl-2p45_kt-1p00_c2-0p00_TuneCP5_13p6TeV_powheg-pythia8"


config.Data.secondaryInputDataset = "/GluGlutoHHto4B_kl-5p00_kt-1p00_c2-0p00_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22EEDR-Poisson70KeepRAW_124X_mcRun3_2022_realistic_postEE_v1-v2/GEN-SIM-RAW"
config.Data.inputDataset          = "/GluGlutoHHto4B_kl-5p00_kt-1p00_c2-0p00_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22EEMiniAODv3-Poisson70KeepRAW_124X_mcRun3_2022_realistic_postEE_v1-v2/MINIAODSIM"
config.Data.outputDatasetTag      = "GluGlutoHHto4B_kl-5p00_kt-1p00_c2-0p00_TuneCP5_13p6TeV_powheg-pythia8"
config.General.requestName        = "GluGlutoHHto4B_kl-5p00_kt-1p00_c2-0p00_TuneCP5_13p6TeV_powheg-pythia8"


config.Data.unitsPerJob  = 10

config.Data.allowNonValidInputDataset = True
config.Data.splitting   = 'EventAwareLumiBased' #Automatic' #'FileBased' #EventAwareLumiBased'
config.Data.publication = False
config.Data.totalUnits  = -1
config.Data.outLFNDirBase = '/store/user/mkolosov/CRAB3_TransferData/HLTRun3/Feb15/'

# Site    
config.Site.storageSite   = 'T3_US_FNALLPC'


