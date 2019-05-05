#!/usr/bin/env python

# Copyright INRA (Institut National de la Recherche Agronomique)
# http://www.inra.fr
# http://urgi.versailles.inra.fr
#
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software.  You can  use,
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info".
#
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability.
#
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or
# data to be ensured and,  more generally, to use and operate it in the
# same conditions as regards security.
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.

import os
import sys
import shutil

sys.path.append(os.environ["REPET_PATH"])
if not "REPET_PATH" in os.environ.keys():
    print "ERROR: no environment variable REPET_PATH"
    sys.exit(1)
sys.path.append(os.environ["REPET_PATH"])
if not "PYTHONPATH" in os.environ.keys():
    os.environ["PYTHONPATH"] = os.environ["REPET_PATH"]
else:
    os.environ["PYTHONPATH"] = "%s:%s" % (os.environ["REPET_PATH"], os.environ["PYTHONPATH"])

from commons.core.LoggerFactory import LoggerFactory
from commons.core.utils.FileUtils import FileUtils
from commons.core.utils.RepetOptionParser import RepetOptionParser
from commons.core.checker.CheckerUtils import CheckerUtils
from commons.core.checker.ConfigChecker import ConfigRules
from commons.core.checker.ConfigChecker import ConfigChecker
from commons.core.checker.CheckerException import CheckerException
from commons.core.seq.FastaUtils import FastaUtils
from commons.core.sql.DbFactory import DbFactory
from commons.core.sql.TableJobAdaptatorFactory import TableJobAdaptatorFactory
from commons.core.launcher.Launcher import Launcher
from commons.tools.LaunchPASTEC import LaunchPASTEC
from commons.tools.GetClassifUniq import GetClassifUniq
from commons.tools.DetectTEFeatures import DetectTEFeatures
from commons.tools.RenameHeaderClassif import RenameHeaderClassif
from commons.tools.NoCatBestHitClassifier import NoCatBestHitClassifier
from commons.tools.RemoveRedundancyBasedOnCI import RemoveRedundancyBasedOnCI
from commons.tools.DetectTEFeatures_parallelized import DetectTEFeatures_parallelized
from commons.tools.ReverseComplementAccordingToClassif import ReverseComplementAccordingToClassif
from PASTEC.StatPastec import StatPastec

LOG_DEPTH = "repet.tools"

####PASTEClassifier
#
class PASTEClassifier(object):

    def __init__(self, fastaFileName = "", configFileName = "", decisionRulesFileName = "", steps = "0", parallel = False, doClean = False, verbosity = 3):
        self._fastaFileName = fastaFileName
        self._addWickerCode = False
        self._reverseComp = False
        self._configFileName = configFileName
        self._iConfig = None

        self._decisionRulesFileName = decisionRulesFileName
        self._steps = steps
        self._parallel = parallel
        self._doClean = doClean
        self._verbosity = verbosity
        self._removeRedundancy = False
        self._addNoCatBestHitClassif = False
        self._classifFileName = ""
        self._projectNameSuffix = ""

        self._projectName = ""
        self._log = LoggerFactory.createLogger("%s.%s" % (LOG_DEPTH, self.__class__.__name__), self._verbosity)

    def setAttributesFromCmdLine(self):
        usage = "PASTEClassifier.py [options]"
        description = "To detect TE features on sequences, classify them and to give some classification statistics.\n"
        description += "Can rename headers with classification info and Wicker's code at the beginning.\n"
        description += "Can reverse-complement consensus if they are detected in reverse strand.\n"
        description += "Warning : it's highly advised to use sequences in upper case.\n"
        epilog = "Example 1: launch in parallel and clean temporary files\n"
        epilog += "\t$ PASTEClassifier.py -i consensus.fa -C PASTEClassifier.cfg -p -c\n"
        epilog += "\n"
        epilog += "Example 2: launch without clean temporary files\n"
        epilog += "\t$ PASTEClassifier.py -i consensus.fa -C PASTEClassifier.cfg \n"
        parser = RepetOptionParser(description = description, epilog = epilog, usage = usage)
        parser.add_option("-i", "--fasta",        dest = "fastaFileName",         action = "store", type = "string",  help = "input fasta file name [compulsory] [format: fasta]", default = "")
        parser.add_option("-C", "--config",       dest = "configFileName",        action = "store", type = "string",  help = "configuration file name (e.g. PASTEClassifier.cfg) [compulsory]", default = "")
        parser.add_option("-D", "--decisionRules",dest = "decisionRulesFileName", action = "store", type = "string",  help = "classification rules file name (e.g. PASTEClassifierRules.yml) [optional]", default = "")
        parser.add_option("-S", "--step",         dest = "step",                  action = "store", type = "string",  help = "step (0/1/2): default: 0 for all steps", default = "0")
        parser.add_option("-p", "--parallel",     dest = "parallel",              action = "store_true",              help = "run tool in parallel", default = False)
        parser.add_option("-c", "--clean",        dest = "doClean",               action = "store_true",              help = "clean temporary files [optional] [default: False]", default = False)
        parser.add_option("-v", "--verbosity",    dest = "verbosity",             action = "store", type = "int",     help = "verbosity [optional] [default: 3, from 1 to 4]", default = 3)
        options = parser.parse_args()[0]
        self._setAttributesFromOptions(options)

    def _setAttributesFromOptions(self, options):
        self.setFastaFileName(options.fastaFileName)
        self.setConfigFileName(options.configFileName)
        self.setDecisionRulesFileName(options.decisionRulesFileName)
        self._steps = options.step
        self._parallel = options.parallel
        self.setDoClean(options.doClean)
        self.setVerbosity(options.verbosity)

    def _checkConfig(self):
        iConfigRules = ConfigRules()
        iConfigRules.addRuleOption(section = "project", option = "project_name", mandatory = True, type = "string")
        sectionName = "classif_consensus"
        iConfigRules.addRuleOption(section = sectionName, option = "clean", mandatory = True, type = "bool")
        iConfigRules.addRuleOption(section = sectionName, option = "remove_redundancy", mandatory = True, type = "bool")
        iConfigRules.addRuleOption(section = sectionName, option = "rev_complement", mandatory = True, type = "bool")
        iConfigRules.addRuleOption(section = sectionName, option = "add_wicker_code", mandatory = True, type = "bool")
        iConfigRules.addRuleOption(section = sectionName, option = "add_noCat_bestHitClassif", mandatory = True, type = "bool")

        if self._parallel:
            iConfigRules.addRuleOption(section = sectionName, option = "limit_job_nb", type = "int")
            iConfigRules.addRuleOption(section = sectionName, option = "resources", type = "string")
            iConfigRules.addRuleOption(section = sectionName, option = "tmpDir", type = "string")

        iConfigChecker = ConfigChecker(self._configFileName, iConfigRules)
        self._iConfig = iConfigChecker.getConfig()
        self._setAttributesFromConfig()

    def _setAttributesFromConfig(self):
        if self._projectName == "":
            self.setProjectName(self._iConfig.get("project", "project_name"))
        projectNameWithoutSuffix = self._projectName.split("_sim")[0]
        projectNameWithoutSuffix = projectNameWithoutSuffix.split("_struct")[0]
        if not CheckerUtils.isMax15Char(projectNameWithoutSuffix):
            print "ERROR: project name '%s' is too long. It must have 15 characters max" % self._projectName
            sys.exit(1)
        if not CheckerUtils.isCharAlphanumOrUnderscore(projectNameWithoutSuffix):
            print "ERROR: project name '%s' must contain only alphanumeric or underscore characters" % self._projectName
            sys.exit(1)


        sectionName = "classif_consensus"
        self.setDoClean(self._iConfig.get(sectionName, "clean"))
        self.setRmRdd(self._iConfig.get(sectionName, "remove_redundancy"))
        self._reverseComp = self._iConfig.get(sectionName, "rev_complement")
        self._addWickerCode = self._iConfig.get(sectionName, "add_wicker_code")
        self._addNoCatBestHitClassif = self._iConfig.get(sectionName, "add_noCat_bestHitClassif")

        if self._parallel:
            self._maxJobNb = self._iConfig.get(sectionName, "limit_job_nb")
            self._resources = self._iConfig.get(sectionName, "resources")
            self._tmpDir = self._iConfig.get(sectionName, "tmpDir")

        self._classifFileName = "%s.classif" % self._projectName

    def setFastaFileName(self, fastaFileName):
        self._fastaFileName = fastaFileName

    def setConfigFileName(self, configFileName):
        self._configFileName = configFileName

    def setDecisionRulesFileName(self, decisionRules):
        self._decisionRulesFileName = decisionRules

    def setAddWickerCode(self, addWickerCode):
        self._addWickerCode = addWickerCode

    def setReverseComp(self, reverseComp):
        self._reverseComp = reverseComp

    def setDoClean(self, doClean):
        self._doClean = doClean

    def setRmRdd(self, removeRedundancy):
        self. _removeRedundancy = removeRedundancy

    def setVerbosity(self, verbosity):
        self._verbosity = verbosity

    def setProjectName(self, projectName):
        self._projectName = projectName

    def _checkOptions(self):
        if self._fastaFileName == "":
            self._logAndRaise("ERROR: Missing input fasta file name")
        else:
            separator = "\n"
            inGenomeFileHandler = open(self._fastaFileName, "r")
            try:
                CheckerUtils.checkHeaders(inGenomeFileHandler)
            except CheckerException, e:
                print "Error in file %s. Wrong headers are :" % self._fastaFileName
                print separator.join(e.messages)
                print "Authorized characters are : a-z A-Z 0-9 - . : _\n"
                inGenomeFileHandler.close()
                sys.exit(1)
            inGenomeFileHandler.close()


    def _logAndRaise(self, errorMsg):
        self._log.error(errorMsg)
        raise Exception(errorMsg)

    def getPASTECcommand(self, iLauncher, fileName):
        lArgs = []
        lArgs.append("-C %s" % self._configFileName)
        if self._decisionRulesFileName:
            lArgs.append("-D %s" % self._decisionRulesFileName)
        lArgs.append("-P %s" % self._projectName)
        lArgs.append("-S 2")
        lArgs.append("-i %s" % fileName)
        lArgs.append("-v %s" % self._verbosity)
        return iLauncher.getSystemCommand("LaunchPASTEC.py", lArgs)

    def _classify(self):
        iLP = LaunchPASTEC(configFileName = self._configFileName, decisionRulesFileName = self._decisionRulesFileName, inputFileName = self._fastaFileName, projectName = self._projectName, verbose = self._verbosity)
        iLP.run()

    def _classifyInParallel(self, nbSeq):
        self._log.debug("Insert banks in database")
        iLP = LaunchPASTEC(configFileName = self._configFileName, step = "1", inputFileName = self._fastaFileName, projectName = self._projectName, verbose = self._verbosity)
        iLP.run()

        self._log.info("Split fasta file")
        minSeqPerJob = 100
        if self._maxJobNb == 0 or nbSeq / self._maxJobNb <= 1.0:
            nbSeqPerBatch = minSeqPerJob
        else:
            nbSeqPerBatch = nbSeq / self._maxJobNb + 1
        FastaUtils.dbSplit(self._fastaFileName, nbSeqPerBatch, True, verbose = self._verbosity - 2)

        self._log.info("Launch PASTEC on each batch")
        queue = self._resources
        cDir = os.getcwd()
        if self._tmpDir != "":
            tmpDir = self._tmpDir
        else:
            tmpDir = cDir

        groupid = "%s_PASTEC" % self._projectName
        acronym = "PASTEC"
        iDb = DbFactory.createInstance()
        iTJA = TableJobAdaptatorFactory.createInstance(iDb, "jobs")
        iLauncher = Launcher(iTJA, cDir, tmpDir, queue, groupid)
        lCmdsTuples = []
        lFiles = FileUtils.getFileNamesList("%s/batches" % cDir, "batch_")
        if len(lFiles) == 0:
            self._logAndRaise("ERROR: directory 'batches' is empty")

        classifFileName = self._classifFileName

        count = 0
        for f in lFiles:
            count += 1
            lCmds = [self.getPASTECcommand(iLauncher, f)]
            lCmdStart = []
            lCmdStart.append("shutil.copy(\"%s/batches/%s\", \".\")" % (cDir, f))
            lCmdStart.append("shutil.copy(\"%s/%s\", \".\")" % (cDir, self._configFileName))
            lCmdFinish = []
            lCmdFinish.append("shutil.move(\"%s\", \"%s/%s_%i\")" % (classifFileName, cDir, classifFileName, count))
            lCmdsTuples.append(iLauncher.prepareCommands_withoutIndentation(lCmds, lCmdStart, lCmdFinish))
        iLauncher.runLauncherForMultipleJobs(acronym, lCmdsTuples, self._doClean)

        FileUtils.catFilesByPattern("%s_*" % classifFileName, classifFileName)
        if self._doClean:
            FileUtils.removeFilesByPattern("%s_*" % classifFileName)
            shutil.rmtree("batches")

    def _postProcessClassification(self):
        self._log.info("Started post processing of Classification")
        newFastaFileName = self._fastaFileName
        newClassifFileName = self._classifFileName

        if os.path.exists(newClassifFileName):
            self._log.debug("Compute stats about classification on initial classif file")
            iSP = StatPastec(inFileName = newClassifFileName)
            iSP.run()

            if self. _removeRedundancy:
                self._log.info("Removing redundancy")
                withoutRddyFastaFileName = "%s_withoutRedundancy.fa" % self._projectName
                iRemoveRedundancy = RemoveRedundancyBasedOnCI(newFastaFileName, newClassifFileName, self._configFileName, outFileName = withoutRddyFastaFileName, doClean = self._doClean, verbosity = self._verbosity)
                iRemoveRedundancy.run()

                # Update classif file after redundancy removal
                iGetClassifuniq = GetClassifUniq(withoutRddyFastaFileName, newClassifFileName, verbosity = self._verbosity)
                iGetClassifuniq.run()
                withoutRddyClassifFileName = "%s.classif" % os.path.splitext(os.path.basename(withoutRddyFastaFileName))[0]

                # Gives stats on classification after redundancy removal
                iSP = StatPastec(inFileName = withoutRddyClassifFileName)
                iSP.run()
                newFastaFileName = withoutRddyFastaFileName
                newClassifFileName = withoutRddyClassifFileName

            if self._reverseComp:
                self._log.info("Reverse complement")
                iRevComplAccording2Classif = ReverseComplementAccordingToClassif(fastaFile = newFastaFileName, classifFileName = newClassifFileName, isOutClassif = True)
                iRevComplAccording2Classif.run()

                newFastaFileName = "%s_negStrandReversed.fa" % os.path.splitext(newFastaFileName)[0]
                newClassifFileName = "%s.classif" % os.path.splitext(newFastaFileName)[0]

            if self._addWickerCode:
                self._log.info("Rename headers according to Wicker's code")
                projectNameInConfigFile = self._iConfig.get("project", "project_name")
                iRHC = RenameHeaderClassif(classifFileName = newClassifFileName, fastaFileName = newFastaFileName, projectName = projectNameInConfigFile, isShorten = True, renameInClassif = True)
                iRHC.run()

                newFastaFileName = "%s_WickerH.fa" % os.path.splitext(newFastaFileName)[0]
                newClassifFileName = "%s.classif" % os.path.splitext(newFastaFileName)[0]

            if self._addNoCatBestHitClassif:
                self._log.info("Adding noCat classification based on bestHit")
                iNCBHC = NoCatBestHitClassifier(fastaFileName = newFastaFileName, classifFileName = newClassifFileName)
                iNCBHC.run()

                newFastaFileName = "%s_noCatBestHit.fa" % os.path.splitext(newFastaFileName)[0]
                newClassifFileName = "%s.classif" % os.path.splitext(newFastaFileName)[0]

            self._log.debug("Compute stats about classification on final classif file (after all post processing)")
            iSP = StatPastec(inFileName = newClassifFileName)
            iSP.run()

            iDb = DbFactory.createInstance()
            iDb.createTable("%s_consensus_classif" % self._projectName, "classif", newClassifFileName, True)
            iDb.close()

#            shutil.move(newFastaFileName, "%s_denovoLibTEs.fa" % self._projectName)
            os.symlink(newFastaFileName, "%s_denovoLibTEs.fa" % self._projectName, )

        else:
            self._logAndRaise("No classification file found or generated")
        self._log.info("Finished post processing of Classification")

    # # Setup the required environment.
    #
    # @param config ConfigParser instance
    #
    def _setup_env(self):
        os.environ["REPET_HOST"] = self._iConfig.get("repet_env", "repet_host")
        os.environ["REPET_USER"] = self._iConfig.get("repet_env", "repet_user")
        os.environ["REPET_PW"] = self._iConfig.get("repet_env", "repet_pw")
        os.environ["REPET_DB"] = self._iConfig.get("repet_env", "repet_db")
        os.environ["REPET_PORT"] = self._iConfig.get("repet_env", "repet_port")
        os.environ["REPET_JOB_MANAGER"] = self._iConfig.get("repet_env", "repet_job_manager")
        os.environ["REPET_QUEUE"] = self._iConfig.get("repet_env", "repet_job_manager")
        os.environ["REPET_JOBS"] = "MySQL"

    def run(self):
        LoggerFactory.setLevel(self._log, self._verbosity)
        if self._configFileName:
            self._checkConfig()
        self._checkOptions()
        self._setup_env()

        toolName = "PASTEClassifier"
        if self._parallel:
            toolName = "PASTEClassifier parallelized"

        self._log.info("START %s" % toolName)
        self._log.info("Fasta file name: %s" % self._fastaFileName)
        nbSeq = FastaUtils.dbSize(self._fastaFileName)
        self._log.debug("Total number of sequences: %i" % nbSeq)

        if "1" in self._steps or "0" in self._steps:
            self._log.info("Running STEP 1 of %s: DetectTEFeatures" % toolName)
            if self._parallel:
                iDF = DetectTEFeatures_parallelized(self._fastaFileName, self._projectName, self._configFileName, self._doClean, self._verbosity)
                iDF.run()
            else:
                iDF = DetectTEFeatures(self._fastaFileName, self._projectName, self._configFileName, self._doClean, self._verbosity)
                iDF.run()
            self._log.info("Finished STEP 1 of %s: DetectTEFeatures" % toolName)

        if "2" in self._steps or "0" in self._steps:
            self._log.info("Running STEP 2 of %s: Classification" % toolName)
            if self._parallel:
                self._classifyInParallel(nbSeq)
            else:
                self._classify()
            self._postProcessClassification()

            self._log.info("Finished STEP 2 of %s: Classification" % toolName)

        self._log.info("END %s" % toolName)

if __name__ == "__main__":
    iPASTEClassifier = PASTEClassifier()
    iPASTEClassifier.setAttributesFromCmdLine()
    iPASTEClassifier.run()
