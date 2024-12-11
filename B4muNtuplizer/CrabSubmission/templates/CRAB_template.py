from CRABClient.UserUtilities import config, getUsername
config = config()

config.General.requestName = 'SkimB4Mu_YEAReraERANAME_Mini'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'

config.JobType.psetName = 'FILE_TO_SUBMIT_PATH'

config.Data.inputDataset = '/DATASET_ID/MINIAOD'
config.Data.allowNonValidInputDataset = True
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 100
config.Data.lumiMask = 'https://cms-service-dqmdc.web.cern.ch/CAF/certification/GOLDEN_JSON_PATH'
#config.Data.publication = True
config.Data.outputDatasetTag = 'SkimB4Mu_YEAReraERANAME_Mini'
config.JobType.allowUndistributedCMSSW = True 
config.Site.storageSite = 'T2_IT_Bari'
config.Site.ignoreGlobalBlacklist  = True
