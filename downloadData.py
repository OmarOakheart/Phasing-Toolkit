import sys
import subprocess
import os
import gzip
import shutil

def updateLog(logFilePath,logText):
    logFile=open(logFilePath,"a")
    logFile.write(logText)
    logFile.close()

def loadBenchmarkParameterFile(parameterFilePath):
    benchmarkParameters = {"mainPath": "",
                           "longReads": {},
                           "longReadMethod": "",
                           "shortReads": {},
                           "coverages": [],
                           "heterozygosityRates": [],
                           "reference": {},
                           "genomeSize": 0,
                           "strainLists": [],
                           "threads":0}
    benchmarkParameterFile=open(parameterFilePath,"r")
    for line in benchmarkParameterFile:
        line=line.strip("\n").split("\t")
        if line[0] in {"mainPath","longReadMethod","genomeSize","threads"}:
            benchmarkParameters[line[0]]=line[1]
        elif line[0] in {"coverages","heterozygosityRates"}:
            for dataPoint in line[1:]:
                benchmarkParameters[line[0]].append(dataPoint)
        elif line[0] in {"longReads","shortReads","reference"}:
            for dataPoint in line[1:]:
                dataPoint=dataPoint.split(":")
                benchmarkParameters[line[0]][dataPoint[0]]=dataPoint[1]
        elif line[0] in {"strainLists"}:
            for dataPoint in line[1:]:
                dataPoint = dataPoint.split(",")
                benchmarkParameters[line[0]].append(dataPoint)
        else:
            print("Ignoring '","\t".join(line),"' because it is not recognized",sep="")
    benchmarkParameterFile.close()
    return benchmarkParameters

def downloadData(benchmarkParameters):
    # Make paths
    mainPath=benchmarkParameters["mainPath"]
    logPath=mainPath+"log/"
    allPaths=[logPath]

    if "shortReads" in benchmarkParameters.keys():
        shortReadPath=mainPath+"rawData/shortReads/"
        allPaths.append(shortReadPath)
    if "longReads" in benchmarkParameters.keys():
        longReadPath=mainPath+"rawData/longReads/"
        allPaths.append(longReadPath)
    if "reference" in benchmarkParameters.keys():
        speciesName=benchmarkParameters["reference"]["speciesName"]
        referencePath=mainPath+"reference/"+speciesName+"/"

    for path in allPaths:
        os.makedirs(path, exist_ok=True)

    # Prepare log
    fullLogPath = logPath + "downloadDataLog.txt"
    logText = ""

    if "shortReads" in benchmarkParameters.keys():
        for strainName, ERRCode in benchmarkParameters["shortReads"].items():
            downloadERR(ERRCode,shortReadPath,fullLogPath)
    if "longReads" in benchmarkParameters.keys():
        for strainName, ERRCode in benchmarkParameters["longReads"].items():
            downloadERR(ERRCode,longReadPath,fullLogPath)
    if "reference" in benchmarkParameters.keys():
        downloadReference(benchmarkParameters,referencePath,fullLogPath)

def downloadERR(ERRCode,dataPath,fullLogPath):
    finalFileNameSR1 = os.path.join(dataPath, ERRCode + "_1.fastq.gz")
    finalFileNameSR2 = os.path.join(dataPath, ERRCode + "_2.fastq.gz")
    finalFileNameLR = os.path.join(dataPath, ERRCode + ".fastq.gz")
    if (os.path.exists(finalFileNameSR1) and os.path.exists(finalFileNameSR2)) or os.path.exists(finalFileNameLR):
        return ERRCode + "already downloaded.\n"
    downloadERRCommand=["sra-downloader",ERRCode,"--save-dir",dataPath]
    p=subprocess.run(downloadERRCommand,stderr=subprocess.PIPE,stdout=subprocess.PIPE, universal_newlines=True)
    logText="COMMAND: "+" ".join(downloadERRCommand)+"\n\n"
    if p.stderr!="" or p.stdout !="":
        logText+="STDERR:\n\n"+p.stderr+"\n\nSTDOUT:\n\n"+p.stdout+"\n\n"
    updateLog(fullLogPath,logText)

def downloadReference(benchmarkParameters,referencePath,fullLogPath):
    taxID = benchmarkParameters["reference"]["taxID"]
    group = benchmarkParameters["reference"]["group"]
    speciesName=benchmarkParameters["reference"]["speciesName"]
    referenceFilePath=referencePath+speciesName+".fasta"

    # Download the reference sequence
    print("Downloading reference sequence")
    p = subprocess.run(["ncbi-genome-download", "-o", referencePath, "--flat-output", "--formats", "fasta", "--taxids", taxID, group], stderr=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
    logText = "COMMAND: " + " ".join(["ncbi-genome-download", "-o", referencePath, "--flat-output", "--formats", "fasta", "--taxids", taxID, group]) + "\n\n"
    if p.stderr != "" or p.stdout != "":
        logText += "STDERR:\n\n" + p.stderr + "\n\nSTDOUT:\n\n" + p.stdout + "\n\n"
    updateLog(fullLogPath, logText)

    gzippedReferenceName = os.listdir(referencePath)[0]

    with gzip.open(referencePath + gzippedReferenceName, 'rb') as f_in:
        with open(referenceFilePath, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

    os.remove(referencePath + gzippedReferenceName)

    # Make sure the reference is indexed
    print("Indexing reference")

    p = subprocess.run(["samtools", "faidx", referenceFilePath], stderr=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
    logText = "COMMAND: " + " ".join(["samtools", "faidx", referenceFilePath]) + "\n\n"
    if p.stderr != "" or p.stdout != "":
        logText += "STDERR:\n\n" + p.stderr + "\n\nSTDOUT:\n\n" + p.stdout + "\n\n"
    updateLog(fullLogPath, logText)

    p = subprocess.run(["bwa", "index", referenceFilePath], stderr=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
    logText = "COMMAND: " + " ".join(["bwa", "index", referenceFilePath]) + "\n\n"
    if p.stderr != "" or p.stdout != "":
        logText += "STDERR:\n\n" + p.stderr + "\n\nSTDOUT:\n\n" + p.stdout + "\n\n"
    updateLog(fullLogPath, logText)

    p = subprocess.run(["gatk", "CreateSequenceDictionary", "-R", referenceFilePath], stderr=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
    logText = "COMMAND: " + " ".join(["gatk", "CreateSequenceDictionary", "-R", referenceFilePath]) + "\n\n"
    if p.stderr != "" or p.stdout != "":
        logText += "STDERR:\n\n" + p.stderr + "\n\nSTDOUT:\n\n" + p.stdout + "\n\n"
    updateLog(fullLogPath, logText)

    print("Reference downloaded and indexed")



if __name__=="__main__":
    # Load parameters
    benchmarkParameterFilePath = sys.argv[1]
    benchmarkParameters = loadBenchmarkParameterFile(benchmarkParameterFilePath)

    # Download data
    downloadData(benchmarkParameters)

    print("Downloads complete")
