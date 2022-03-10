import sys
import os
import subprocess
import timeit
import psutil
import time

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

def launchFunction(functionName,parameters):
    if functionName=="nphase":
        strainName=parameters[0]
        referenceFilePath=parameters[1]
        outputFolder=parameters[2]
        longReadFile=parameters[3]
        vcfFile=parameters[4]
        mappedLR=parameters[5]
        threadNumber=parameters[6]
        outputLog, systemMessage=launchnPhase(strainName, referenceFilePath, outputFolder, longReadFile, vcfFile, mappedLR, threadNumber)
        print(systemMessage)
        updateLog(fullLogPath, outputLog)
    elif functionName=="whatshap Polyphase":
        ploidy=parameters[0]
        whatsHapPolyphaseParameter=parameters[1]
        outPath=parameters[2]
        referenceFilePath=parameters[3]
        threads=parameters[4]
        vcfFile=parameters[5]
        mappedLR=parameters[6]
        indelBool=parameters[7]
        localStrainName=parameters[8]
        outputLog, systemMessage=launchWhatshapPolyphase(ploidy,whatsHapPolyphaseParameter,outPath,referenceFilePath,threads,
                                                         vcfFile,mappedLR,indelBool,localStrainName)
        print(systemMessage)
        updateLog(fullLogPath, outputLog)
    elif functionName=="flopp":
        mappedLR=parameters[0]
        vcfFile=parameters[1]
        ploidy=parameters[2]
        outPath=parameters[3]
        threads=parameters[4]
        localStrainName=parameters[5]
        outputLog, systemMessage=launchFlopp(mappedLR,vcfFile,ploidy,outPath,threads,localStrainName)
        print(systemMessage)
        updateLog(fullLogPath, outputLog)

def launchFlopp(mappedLR,vcfFile,ploidy,outPath,threads,localStrainName):
    start = timeit.default_timer()
    outputLog = ""

    p = subprocess.Popen(["flopp","-b",mappedLR, "-c", vcfFile, "-p", ploidy, "-o", outPath, "-t", threads], stderr=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
    outputLog += "COMMAND: " + " ".join(["flopp","-b",mappedLR, "-c", vcfFile, "-p", ploidy, "-o", outPath, "-t", threads]) + "\n\n"

    maxRSS_mem = 0
    maxVMS_mem = 0

    try:
        poll = p.poll()
    except:
        poll = "Done"
    while poll is None:
        monitorP = psutil.Process(p.pid)
        children = list(monitorP.children(recursive=True))
        children = children + [monitorP]
        RSS_mem = 0
        VMS_mem = 0
        for child in children:
            mem_info = child.memory_info()
            RSS_mem += mem_info[0]
            VMS_mem += mem_info[1]
            if RSS_mem > maxRSS_mem:
                maxRSS_mem = RSS_mem
            if VMS_mem > maxVMS_mem:
                maxVMS_mem = VMS_mem
        time.sleep(1)
        try:
            poll = p.poll()
        except:
            poll = "Done"

    if p.stderr != "" or p.stdout != "":
        outputLog += "STDERR:\n\n" + p.stderr.read() + "\n\nSTDOUT:\n\n" + p.stdout.read() + "\n\n"

    print("Max memory usage, RSS:", maxRSS_mem / 1024 / 1024 / 1024, "VMS:", maxVMS_mem / 1024 / 1024 / 1024)

    stop = timeit.default_timer()
    totalRunTime = stop - start
    print("This took", totalRunTime, "seconds to run.")
    performanceLine = "\t".join([str(x) for x in [totalRunTime, maxRSS_mem, maxVMS_mem, "flopp", localStrainName]]) + "\n"
    pFile = open(performanceMetricFile, "a")
    pFile.write(performanceLine)
    pFile.close()
    return outputLog, "flopp ran successfully"

def launchnPhase(strainName,referenceFilePath,outputFolder,longReadFile,vcfFile,mappedLR,threadNumber):
    start = timeit.default_timer()
    outputLog = ""

    p = subprocess.Popen(["nphase", "partial", "--sampleName", strainName, "--reference", referenceFilePath, "--output",
                        outputFolder, "--longReads", longReadFile, "--vcf", vcfFile, "--mappedLongReads",
                        mappedLR, "--threads", threadNumber], stderr=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
    outputLog+="COMMAND: "+" ".join(["nphase partial --sampleName",strainName,"--reference",referenceFilePath,"--output",
                                     outputFolder,"--longReads",longReadFile,"--vcf",vcfFile,"--mappedLongReads",
                                     mappedLR,"--threads", threadNumber]) + "\n\n"

    maxRSS_mem=0
    maxVMS_mem=0

    try:
        poll = p.poll()
    except:
        poll="Done"
    while poll is None:
        monitorP=psutil.Process(p.pid)
        children=list(monitorP.children(recursive=True))
        children=children+[monitorP]
        RSS_mem=0
        VMS_mem=0
        for child in children:
            mem_info=child.memory_info()
            RSS_mem+=mem_info[0]
            VMS_mem+=mem_info[1]
            if RSS_mem>maxRSS_mem:
                maxRSS_mem=RSS_mem
            if VMS_mem>maxVMS_mem:
                maxVMS_mem=VMS_mem
        time.sleep(1)
        try:
            poll=p.poll()
        except:
            poll="Done"


    if p.stderr!="" or p.stdout!="":
        outputLog+="STDERR:\n\n"+p.stderr.read()+"\n\nSTDOUT:\n\n"+p.stdout.read()+"\n\n"

    print("Max memory usage, RSS:",maxRSS_mem/1024/1024/1024,"VMS:",maxVMS_mem/1024/1024/1024)

    stop = timeit.default_timer()
    totalRunTime = stop - start
    print("This took", totalRunTime, "seconds to run.")
    performanceLine="\t".join([str(x) for x in [totalRunTime,maxRSS_mem,maxVMS_mem,"nphase",strainName]])+"\n"
    pFile=open(performanceMetricFile,"a")
    pFile.write(performanceLine)
    pFile.close()
    return outputLog, "nPhase ran successfully"

def launchWhatshapPolyphase(ploidy,whatsHapPolyphaseParameter,outPath,referenceFilePath,threads,vcfFile,mappedLR,indelBool,localStrainName):
    start = timeit.default_timer()

    outputLog = ""

    if indelBool==False:
        p = subprocess.Popen(["whatshap","polyphase","--ploidy",ploidy,"--ignore-read-groups","--block-cut-sensitivity",
                              whatsHapPolyphaseParameter ,"-o",outPath,"--reference",referenceFilePath,"--include-haploid-sets","--threads",threads,
                              vcfFile,mappedLR], stderr=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
        outputLog += "COMMAND: " + " ".join(["whatshap","polyphase","--ploidy",ploidy,"--ignore-read-groups","--block-cut-sensitivity",
                                             whatsHapPolyphaseParameter ,"-o",outPath,"--reference",referenceFilePath,"--include-haploid-sets","--threads",
                                             threads,vcfFile,mappedLR]) + "\n\n"
    else:
        p = subprocess.Popen(["whatshap", "polyphase", "--ploidy", ploidy, "--indels","--ignore-read-groups",
                              "--block-cut-sensitivity",whatsHapPolyphaseParameter, "-o", outPath, "--reference",
                              referenceFilePath,"--include-haploid-sets", "--threads", threads,vcfFile, mappedLR], stderr=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
        outputLog += "COMMAND: " + " ".join(["whatshap", "polyphase", "--ploidy", ploidy,"--indels", "--ignore-read-groups",
                                             "--block-cut-sensitivity",whatsHapPolyphaseParameter, "-o", outPath, "--reference",
                                             referenceFilePath,"--include-haploid-sets", "--threads",threads, vcfFile, mappedLR]) + "\n\n"

    maxRSS_mem = 0
    maxVMS_mem = 0

    try:
        poll = p.poll()
    except:
        poll = "Done"
    while poll is None:
        monitorP = psutil.Process(p.pid)
        children = list(monitorP.children(recursive=True))
        children = children + [monitorP]
        RSS_mem = 0
        VMS_mem = 0
        for child in children:
            mem_info = child.memory_info()
            RSS_mem += mem_info[0]
            VMS_mem += mem_info[1]
            if RSS_mem > maxRSS_mem:
                maxRSS_mem = RSS_mem
            if VMS_mem > maxVMS_mem:
                maxVMS_mem = VMS_mem
        time.sleep(1)
        try:
            poll = p.poll()
        except:
            poll = "Done"

    if p.stderr != "" or p.stdout != "":
        outputLog += "STDERR:\n\n" + p.stderr.read() + "\n\nSTDOUT:\n\n" + p.stdout.read() + "\n\n"

    print("Max memory usage, RSS:", maxRSS_mem / 1024 / 1024 / 1024, "VMS:", maxVMS_mem / 1024 / 1024 / 1024)

    stop = timeit.default_timer()
    totalRunTime = stop - start
    print("This took", totalRunTime / 60 / 60, "hours to run.")
    performanceLine="\t".join([str(x) for x in [totalRunTime,maxRSS_mem,maxVMS_mem,"whatshap_polyphase",localStrainName]])+"\n"
    pFile=open(performanceMetricFile,"a")
    pFile.write(performanceLine)
    pFile.close()
    return outputLog, "whatshap polyphase ran successfully"

if __name__ == "__main__":
    # Load parameters
    benchmarkParameterFilePath = sys.argv[1]
    benchmarkParameters = loadBenchmarkParameterFile(benchmarkParameterFilePath)
    phasingMethods=sys.argv[2:]

    # Handle paths
    mainPath = benchmarkParameters["mainPath"]

    speciesName=benchmarkParameters["reference"]["speciesName"]

    predictionPath=mainPath+"phasingPredictions/"
    hybridPath=mainPath+"virtualPolyploids/"
    hybridRawReads=hybridPath+"rawReads/"
    hybridLongReads=hybridRawReads+"longReads/"
    referencePath=mainPath+"reference/"+speciesName+"/"
    referenceFilePath=referencePath+speciesName+".fasta"
    hybridPath=mainPath+"virtualPolyploids/"
    mappedLongReads=hybridPath+"Mapped/longReads/"
    variantCalledShortReads=hybridPath+"Variants/shortReads/"
    performanceMetricPath=mainPath+"performanceMetrics/"

    allPaths=[predictionPath,performanceMetricPath]

    if "flopp" in phasingMethods:
        floppPath = predictionPath + "flopp/"
        allPaths.append(floppPath)
    if "nphase" in phasingMethods:
        nPhasePath=predictionPath+"nPhase/"
        allPaths.append(nPhasePath)
    if "whatshap-polyphase" in phasingMethods:
        whatsHapPolyphasePath=predictionPath+"whatsHapPolyphase/"
        allPaths.append(whatsHapPolyphasePath)

    logPath=mainPath+"log/"

    for path in allPaths:
        os.makedirs(path, exist_ok=True)

    #Prepare log

    fullLogPath=logPath+"phasingLog.txt"
    logText=""

    #Prepare performance metrics file
    performanceMetricFile=performanceMetricPath+"timeAndMemoryMetrics.txt"
    pFile = open(performanceMetricFile, "a")
    pFile.write("#time (seconds)\tRSS (bytes)\tVMS (bytes)\ttoolName\ttestName\n")
    pFile.close()

    threads=benchmarkParameters["threads"]

    #Run all flopp
    if "flopp" in phasingMethods:
        for strainList in benchmarkParameters["strainLists"]:
            strainName="_".join(strainList)
            ploidy = str(len(strainList))
            for coverage in benchmarkParameters["coverages"]:
                for heterozygosityRate in benchmarkParameters["heterozygosityRates"]:
                    testName = strainName + "_" + str(coverage) + "X_"+heterozygosityRate
                    vcfFile=variantCalledShortReads+testName+".SNPs.vcf"
                    mappedLR = mappedLongReads + strainName + "_" + str(coverage) + "X.sorted.bam"
                    outPath = floppPath+"flopp_"+testName+".SNPs.vcf"
                    launchFunction("flopp",[mappedLR, vcfFile, ploidy, outPath, threads,testName])
                    vcfFileWithIndels = variantCalledShortReads + testName + ".vcf"
                    outPath = floppPath + "flopp_" + testName + ".vcf"
                    launchFunction("flopp", [mappedLR, vcfFileWithIndels, ploidy, outPath, threads,testName+"_Indels"])

        print("Done with all flopp")

    #Run all whatshap polyphase
    if "whatshap-polyphase" in phasingMethods:
        for strainList in benchmarkParameters["strainLists"]:
            strainName="_".join(strainList)
            ploidy = str(len(strainList))
            for coverage in benchmarkParameters["coverages"]:
                for heterozygosityRate in benchmarkParameters["heterozygosityRates"]:
                    testName = strainName + "_" + str(coverage) + "X_" + heterozygosityRate
                    vcfFile=variantCalledShortReads+testName+".SNPs.vcf"
                    mappedLR = mappedLongReads + strainName + "_" + str(coverage) + "X.sorted.bam"
                    outPath = whatsHapPolyphasePath+"WHP_"+testName+".SNPs.vcf"
                    whatsHapPolyphaseParameter = "4"
                    indelBool = False
                    launchFunction("whatshap Polyphase",[ploidy,whatsHapPolyphaseParameter,outPath,referenceFilePath,threads,vcfFile,
                                                         mappedLR,indelBool,testName])
                    vcfFileWithIndels = variantCalledShortReads + testName + ".vcf"
                    indelBool=True
                    outPath = whatsHapPolyphasePath + "WHP_" + testName + ".vcf"
                    launchFunction("whatshap Polyphase", [ploidy, whatsHapPolyphaseParameter, outPath, referenceFilePath, threads,
                                                          vcfFileWithIndels,mappedLR, indelBool,testName+"_Indels"])

        print("Done with all whatshap polyphase")

    #Run all nPhase
    if "nphase" in phasingMethods:
        for strainList in benchmarkParameters["strainLists"]:
            strainName="_".join(strainList)
            ploidy = str(len(strainList))
            for coverage in benchmarkParameters["coverages"]:
                for heterozygosityRate in benchmarkParameters["heterozygosityRates"]:
                    testName=strainName+"_"+str(coverage)+"X_"+heterozygosityRate
                    longReadFile=hybridLongReads+strainName+"_"+str(coverage)+"X.fastq.gz"
                    vcfFile=variantCalledShortReads+testName+".SNPs.vcf"
                    mappedLR=mappedLongReads+strainName+"_"+str(coverage)+"X.sorted.sam"
                    launchFunction("nphase", [testName,referenceFilePath,nPhasePath,longReadFile,vcfFile,mappedLR,threads])
                    vcfFileWithIndels=variantCalledShortReads+testName+".vcf"
                    launchFunction("nphase", [testName+"_Indels",referenceFilePath,nPhasePath,longReadFile,vcfFileWithIndels,mappedLR,threads])

        print("Done with all nPhase")

    print("All datasets phased with all tools")
