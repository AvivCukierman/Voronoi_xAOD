import ROOT
import logging
import shutil
import os

logging.basicConfig(level=logging.INFO)
from optparse import OptionParser

parser = OptionParser()

# job configuration
parser.add_option("--submitDir", help="dir to store the output", default="submit_dir_z1and2")
parser.add_option("--inputDQ2", dest='use_scanDQ2', action='store_true', help="input files from DQ2", default=False )
parser.add_option("--inputFiles", type=str, dest='search_directories', help="list of input files", default = "/atlas/local/acukierm/dijetz1and2/")
parser.add_option("--driver", help="select where to run", choices=("direct", "prooflite", "grid", "lsf"), default="direct")
parser.add_option("--nevents", type=int, help="number of events to process for all the datasets")
parser.add_option("--skip-events", type=int, help="skip the first n events")
parser.add_option("-w", "--overwrite", action='store_true', default=False, help="overwrite previous submitDir")
# input configuration
parser.add_option("--dataSource", help="undefined=-1, data=0, FullSim=1, AF-II=2 ", type="int", default=1)

(options, args) = parser.parse_args()

import atexit
@atexit.register
def quite_exit():
    ROOT.gSystem.Exit(0)

logging.info("loading packages")
ROOT.gROOT.Macro("$ROOTCOREDIR/scripts/load_packages.C")

if options.overwrite:
    import shutil
    shutil.rmtree(options.submitDir, True)

#Set up the job for xAOD access:
ROOT.xAOD.Init().ignore();

use_scanDQ2 = options.use_scanDQ2 or options.driver in ['grid','lsf']
if use_scanDQ2:
  if os.getenv('XRDSYS') is None:
    raise EnvironmentError('xrootd client is not setup. Run localSetupFAX or equivalent.')

if options.driver == 'lsf':
  if getattr(ROOT.EL, 'LSFDriver') is None:
    raise KeyError('Cannot load the LSF driver from EventLoop. Did you not compile it?')

# create a new sample handler to describe the data files we use
logging.info("creating new sample handler")
sh_all = ROOT.SH.SampleHandler()

#search_directories = ["/afs/cern.ch/work/a/amarzin/susy_13TeV/samples/W/"]
#search_directories = ["/nfs/slac/g/atlas/u01/users/acukierm/xAOD2/MyPackage/scripts"]
'''if use_scanDQ2:
  search_directories = ["mc15_13TeV.361020.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ0W.merge.AOD.e3569_s2576_s2132_r6765_r6282/"]
  #search_directories = ["user.amarzin.407012.ttbar.DAOD_SUSY10.e4023_s2608_r6765_r6282_p2375_tag_06_v3_output_xAOD.root/"]
else:
  search_directories = ["/atlas/local/acukierm/dijetz1and2/"]
'''
# scan for datasets in the given directories
directory = options.search_directories
if use_scanDQ2:
  ROOT.SH.scanDQ2(sh_all, directory)
else:
  ROOT.SH.scanDir(sh_all, directory)

# print out the samples we found
logging.info("%d different datasets found scanning all directories", len(sh_all))

# set the name of the tree in our files
sh_all.setMetaString("nc_tree", "CollectionTree")

# this is the basic description of our job
logging.info("creating new job")
job = ROOT.EL.Job()
job.sampleHandler(sh_all)

if options.driver == 'lsf':
  job.options().setBool(ROOT.EL.Job.optResetShell, False);

#Set the xAOD access mode of the job:
job.options().setString( ROOT.EL.Job.optXaodAccessMode,ROOT.EL.Job.optXaodAccessMode_branch );
        
if options.nevents:
    logging.info("processing only %d events", options.nevents)
    job.options().setDouble(ROOT.EL.Job.optMaxEvents, options.nevents)

if options.skip_events:
    logging.info("skipping first %d events", options.skip_events)
    job.options().setDouble(ROOT.EL.Job.optSkipEvents, options.skip_events)

#define an output and an ntuple associated to that output
'''output = ROOT.EL.OutputStream("myOutput");
job.outputAdd(output)
ntuple = ROOT.EL.NTupleSvc("myOutput");
job.algsAdd(ntuple)'''

# add our algorithm to the job
logging.info("creating algorithms")
alg1 = ROOT.VoronoiWeights()
alg2 = ROOT.VoronoiJets()
alg3 = ROOT.JetMatching()
alg4 = ROOT.WriteTree()

logging.info("adding algorithms")
job.algsAdd(alg1)
job.algsAdd(alg2)
job.algsAdd(alg3)
job.algsAdd(alg4)

# make the driver we want to use:
# this one works by running the algorithm directly
logging.info("creating driver")
driver = None
if (options.driver == "direct"):
  logging.info("running on direct")
  driver = ROOT.EL.DirectDriver()
  logging.info("submit job")
  driver.submit(job, options.submitDir)
elif (options.driver == "prooflite"):
  logging.info("running on prooflite")
  driver = ROOT.EL.ProofDriver()
  logging.info("submit job")
  driver.submit(job, options.submitDir)
elif (options.driver == "grid"):
  logging.info("running on Grid")
  driver = ROOT.EL.PrunDriver()   
  driver.options().setString("nc_outputSampleName", "user.%nickname%.%in:name[2]%.%in:name[3]%.%in:name[6]%_11") 
  driver.options().setDouble("nc_nGBPerJob", 2)
  #driver.options().setDouble("nc_disableAutoRetry", 1)
  #driver.options().setString("nc_cmtConfig", "x86_64-slc6-gcc48-opt")
  #driver.options().setDouble("nc_mergeOutput", 1);  
  #driver.options().setString("nc_excludedSite", "TRIUMF-LCG2_DATADISK")
  #driver.options().setString("nc_EventLoop_SubmitFlags", "--allowTaskDuplication");
  driver.options().setString("EL::Job::optGridDestSE","CERN-PROD_LOCALGROUPDISK")
  driver.options().setString("nc_optGridDestSE","CERN-PROD_LOCALGROUPDISK")
  logging.info("submit job")
  driver.submitOnly(job, options.submitDir)
elif (options.driver == "lsf"):
  logging.info("running on LSF")
  driver = ROOT.EL.LSFDriver()
  driver.options().setString(ROOT.EL.Job.optSubmitFlags, '-q short')
  driver.submitOnly(job, options.submitDir)
