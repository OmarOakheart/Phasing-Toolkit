import sys
import os
import subprocess
import random
from multiprocessing import Pool

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

def longReadMapping(strainName,longReads,reference,outputFolder,flag,longReadPlatform,threads):
    outputLog=""
    finalFileName=os.path.join(outputFolder,strainName+".sorted.sam")
    if os.path.exists(finalFileName):
        return outputLog,finalFileName+" already exists.\n"
    #Map with ngmlr
    samFile=os.path.join(outputFolder,strainName+".sam")
    p=subprocess.run(["ngmlr","-t",threads,"-x",longReadPlatform,"-r",reference,"-q",longReads,"-o",samFile],stderr=subprocess.PIPE,stdout=subprocess.PIPE, universal_newlines=True)
    outputLog+="COMMAND: "+" ".join(["ngmlr","-t",threads,"-x",longReadPlatform,"-r",reference,"-q",longReads,"-o",samFile])+"\n\n"
    if p.stderr!="" or p.stdout !="":
        outputLog+="STDERR:\n\n"+p.stderr+"\n\nSTDOUT:\n\n"+p.stdout+"\n\n"

    #Remove unmapped reads
    passSamFile=os.path.join(outputFolder,strainName+".pass.sam")
    p=subprocess.run(["samtools","view","-h","-t",reference,"-@",threads,"-F",flag,samFile,"-o",passSamFile],stderr=subprocess.PIPE,stdout=subprocess.PIPE, universal_newlines=True)
    outputLog+="COMMAND: "+" ".join(["samtools","view","-t",reference,"-@",threads,"-F",flag,samFile,"-o",passSamFile])+"\n\n"
    if p.stderr!="" or p.stdout !="":
        outputLog+="STDERR:\n\n"+p.stderr+"\n\nSTDOUT:\n\n"+p.stdout+"\n\n"

    #Sort sam file
    sortedHeaderSamFile=os.path.join(outputFolder,strainName+".sorted.header.sam")
    p=subprocess.run(["samtools","sort",passSamFile,"-@",threads,"-o",sortedHeaderSamFile],stderr=subprocess.PIPE,stdout=subprocess.PIPE, universal_newlines=True)
    outputLog+="COMMAND: "+" ".join(["samtools","sort",passSamFile,"-@",threads,"-o",sortedHeaderSamFile])+"\n\n"
    if p.stderr!="" or p.stdout !="":
        outputLog+="STDERR:\n\n"+p.stderr+"\n\nSTDOUT:\n\n"+p.stdout+"\n\n"

    #Adding read groups
    sortedHeaderRGSamFile=os.path.join(outputFolder, strainName + "sorted.header.RG.sam")
    p=subprocess.run(["gatk", "AddOrReplaceReadGroups", "-I", sortedHeaderSamFile, "-O", sortedHeaderRGSamFile, "-RGID", "ID_" + strainName, "-RGLB", "LB_" + strainName, "-RGPL", "LONGREADS", "-RGPU", "PU_" + strainName, "-RGSM", "SM_" + strainName], stderr=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
    outputLog+="COMMAND: "+" ".join(["gatk", "AddOrReplaceReadGroups", "-I", sortedHeaderSamFile, "-O", sortedHeaderRGSamFile, "-RGID", "ID_" + strainName, "-RGLB", "LB_" + strainName, "-RGPL", "LONGREADS", "-RGPU", "PU_" + strainName, "-RGSM", "SM_" + strainName]) + "\n\n"
    if p.stderr!=""or p.stdout!="":
        outputLog+="STDERR:\n\n"+p.stderr+"\n\nSTDOUT:\n\n"+p.stdout+"\n\n"

    #Generate bam file
    sortedBamFile = os.path.join(outputFolder, strainName + ".sorted.bam")
    p = subprocess.run(["samtools", "view", "-S", "-b", sortedHeaderRGSamFile, "-@", threads, "-o", sortedBamFile], stderr=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
    outputLog += "COMMAND: " + " ".join(["samtools", "view", "-S", "-b", sortedHeaderRGSamFile, "-@", threads, "-o", sortedBamFile]) + "\n\n"
    if p.stderr != "" or p.stdout != "":
        outputLog += "STDERR:\n\n" + p.stderr + "\n\nSTDOUT:\n\n" + p.stdout + "\n\n"

    #Index bam file
    p = subprocess.run(["samtools", "index", sortedBamFile], stderr=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
    outputLog += "COMMAND: " + " ".join(["samtools", "index",sortedBamFile]) + "\n\n"
    if p.stderr != "" or p.stdout != "":
        outputLog += "STDERR:\n\n" + p.stderr + "\n\nSTDOUT:\n\n" + p.stdout + "\n\n"

    #Get rid of the SAM header (it gets in the way of another function later)
    sortedSamFile=os.path.join(outputFolder,strainName+".sorted.sam")
    p=subprocess.run(["samtools","view",sortedHeaderSamFile,"-@",threads,"-o",sortedSamFile],stderr=subprocess.PIPE,stdout=subprocess.PIPE, universal_newlines=True)
    outputLog+="COMMAND: "+" ".join(["samtools","view",sortedHeaderSamFile,"-@",threads,"-o",sortedSamFile])+"\n\n"
    if p.stderr!="" or p.stdout !="":
        outputLog+="STDERR:\n\n"+p.stderr+"\n\nSTDOUT:\n\n"+p.stdout+"\n\n"

    try:
        #Cleaning up
        os.remove(samFile)
        os.remove(passSamFile)
        os.remove(sortedHeaderSamFile)
        os.remove(sortedHeaderRGSamFile)
    except:
        return outputLog,"Long reads not mapped to reference"

    return outputLog,"Long reads successfully mapped to reference"

def shortReadMapping(strainName,R1,R2,reference,outputFolder,threads):
    outputLog=""
    finalFileName = os.path.join(outputFolder, strainName + ".final.bam")
    if os.path.exists(finalFileName):
        return outputLog, finalFileName + " already exists.\n"
    #Mapping
    samFileNoRGPath=os.path.join(outputFolder, strainName+".noRG.sam")
    samFileNoRG=open(samFileNoRGPath,"w")
    p=subprocess.run(["bwa","mem","-t",threads,"-M", reference, R1, R2],stderr=subprocess.PIPE,stdout=samFileNoRG, universal_newlines=True)
    samFileNoRG.close()
    outputLog+="COMMAND: "+" ".join(["bwa","mem","-t",threads,"-M", reference, R1, R2, ">", samFileNoRGPath])+"\n\n"
    if p.stderr!="":
        outputLog+="STDERR:\n\n"+p.stderr+"\n\nSTDOUT is in "+samFileNoRGPath+"\n\n"

    #Adding read groups
    samFile=os.path.join(outputFolder, strainName+".sam")
    p=subprocess.run(["gatk", "AddOrReplaceReadGroups","-I",samFileNoRGPath,"-O",samFile,"-RGID","ID_"+strainName, "-RGLB","LB_"+strainName, "-RGPL","SHORTREADS", "-RGPU","PU_"+strainName, "-RGSM","SM_"+strainName],stderr=subprocess.PIPE,stdout=subprocess.PIPE, universal_newlines=True)
    outputLog+="COMMAND: "+" ".join(["gatk", "AddOrReplaceReadGroups","-I",samFileNoRGPath,"-O",samFile,"-RGID","ID_"+strainName, "-RGLB","LB_"+strainName, "-RGPL","SHORTREADS", "-RGPU","PU_"+strainName, "-RGSM","SM_"+strainName])+"\n\n"
    if p.stderr!="" or p.stdout!="":
        outputLog+="STDERR:\n\n"+p.stderr+"\n\nSTDOUT:\n\n"+p.stdout+"\n\n"

    #Organizing mapped reads
    bamFile=os.path.join(outputFolder,strainName+".bam")
    sortedBamFile=os.path.join(outputFolder,strainName+".sorted.bam")

    print("Removing quality 0 (multimapped) reads, turning to bam and sorting it")

    p=subprocess.run(["samtools","view","-@",threads,"-bT",reference,"-q","1",samFile,"-o",bamFile],stderr=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
    outputLog+="COMMAND: "+" ".join(["samtools","view","-@",threads,"-bT",reference,"-q","1",samFile,"-o",bamFile])+"\n\n"
    if p.stderr!="" or p.stdout !="":
        outputLog+="STDERR:\n\n"+p.stderr+"\n\nSTDOUT:\n\n"+p.stdout+"\n\n"

    p=subprocess.run(["samtools", "sort", bamFile,"-@",threads, "-o", sortedBamFile],stderr=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
    outputLog+="COMMAND: "+" ".join(["samtools", "sort", bamFile,"-@",threads, "-o", sortedBamFile])+"\n\n"
    if p.stderr!="" or p.stdout !="":
        outputLog+="STDERR:\n\n"+p.stderr+"\n\nSTDOUT:\n\n"+p.stdout+"\n\n"

    p=subprocess.run(["samtools", "index","-@",threads, sortedBamFile],stderr=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
    outputLog+="COMMAND: "+" ".join(["samtools", "index","-@",threads, sortedBamFile])+"\n\n"
    if p.stderr!="" or p.stdout !="":
        outputLog+="STDERR:\n\n"+p.stderr+"\n\nSTDOUT:\n\n"+p.stdout+"\n\n"

    #GATK cleaning
    MDsortedBamFile=os.path.join(outputFolder,strainName+".sorted.MD.bam")
    p=subprocess.run(["gatk","MarkDuplicates", "--REMOVE_DUPLICATES", "true", "-I",sortedBamFile,"-O",MDsortedBamFile,"-M", "/dev/null"],stderr=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
    outputLog+="COMMAND: "+" ".join(["gatk","MarkDuplicates", "--REMOVE_DUPLICATES", "true", "-I",sortedBamFile,"-O",MDsortedBamFile,"-M", "/dev/null"])+"\n\n"
    if p.stderr!="" or p.stdout !="":
        outputLog+="STDERR:\n\n"+p.stderr+"\n\nSTDOUT:\n\n"+p.stdout+"\n\n"

    #Finalizing
    finalBamFile=os.path.join(outputFolder,strainName+".final.bam")
    p=subprocess.run(["samtools", "sort", MDsortedBamFile, "-@",threads,"-o", finalBamFile],stderr=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
    outputLog+="COMMAND: "+" ".join(["samtools", "sort", MDsortedBamFile, "-@",threads,"-o", finalBamFile])+"\n\n"
    if p.stderr!="" or p.stdout !="":
        outputLog+="STDERR:\n\n"+p.stderr+"\n\nSTDOUT:\n\n"+p.stdout+"\n\n"

    p=subprocess.run(["samtools", "index","-@",threads, finalBamFile],stderr=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
    outputLog+="COMMAND: "+" ".join(["samtools", "index","-@",threads, finalBamFile])+"\n\n"
    if p.stderr!="" or p.stdout !="":
        outputLog+="STDERR:\n\n"+p.stderr+"\n\nSTDOUT:\n\n"+p.stdout+"\n\n"

    try:
        #Cleaning up
        os.remove(samFileNoRGPath)
        os.remove(samFile)
        os.remove(bamFile)
        os.remove(sortedBamFile)
        os.remove(MDsortedBamFile)
    except:
        return outputLog,"ERROR: Short reads unsuccessfully mapped to reference"

    return outputLog,"Short reads successfully mapped to reference"

def runVariantCalling(strainName,mappingPath,variantCallingPath,referenceFilePath,estimatedPloidy):
    shortReadBam = os.path.join(mappingPath, strainName + ".final.bam")
    shortReadVCF = os.path.join(variantCallingPath, strainName + ".vcf")
    if os.path.exists(shortReadVCF)==False:
        p = subprocess.run(["gatk", "HaplotypeCaller", "-R", referenceFilePath, "-ploidy", estimatedPloidy, "-I", shortReadBam, "-O", shortReadVCF], stderr=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
        logText = "COMMAND: " + " ".join(["gatk", "HaplotypeCaller", "-R", referenceFilePath, "-ploidy", estimatedPloidy, "-I", shortReadBam, "-O", shortReadVCF]) + "\n\n"
        if p.stderr != "" or p.stdout != "":
            logText += "STDERR:\n\n" + p.stderr + "\n\nSTDOUT:\n\n" + p.stdout + "\n\n"
    else:
        logText=shortReadVCF+" already exists."
    updateLog(fullLogPath, logText)

def subsampleVCF(VCFFile,subsamplingSize,noIndels,referenceFilePath,genomeSize):
    #Only keep N random het SNPs (coverage based). Can be extended to indels later.
    nHetAllowed=int(float(subsamplingSize)/100*int(genomeSize))

    if noIndels:
        SNPsOnlyVCF=VCFFile.replace(".vcf",".SNPs.vcf")
        if os.path.exists(SNPsOnlyVCF)==False:
            p=subprocess.run(["gatk","SelectVariants","-R",referenceFilePath,"--variant",VCFFile,"-O",SNPsOnlyVCF,"--select-type-to-include","SNP"],stderr=subprocess.PIPE,stdout=subprocess.PIPE, universal_newlines=True)
            logText="COMMAND: "+" ".join(["gatk","SelectVariants","-R",referenceFilePath,"--variant",VCFFile,"-O",SNPsOnlyVCF,"--select-type-to-include","SNP"])+"\n\n"
            if p.stderr!="" or p.stdout !="":
                logText+="STDERR:\n\n"+p.stderr+"\n\nSTDOUT:\n\n"+p.stdout+"\n\n"
            updateLog(fullLogPath,logText)
        subsampledVCFFilePath=VCFFile.replace(".vcf","_"+str(subsamplingSize)+".SNPs.vcf")
        SNPVCFFile = open(SNPsOnlyVCF, "r")
    else:
        subsampledVCFFilePath=VCFFile.replace(".vcf","_"+str(subsamplingSize)+".vcf")
        SNPVCFFile = open(VCFFile, "r")

    if os.path.exists(subsampledVCFFilePath):
        print(subsampledVCFFilePath,"already exists.")
    else:
        hetSNPVCFText=""

        allLines=[]
        allHetIndexes=[]
        index=0
        for line in SNPVCFFile:
            if "#" in line:
                allLines.append(line)
            else:
                allLines.append(line)
                format=line.split("\t")[8].split(":")
                if "GT" in format:
                    GTIndex=format.index("GT")
                GTData=line.split("\t")[9].split(":")[GTIndex]
                nHaplotypes=len(set(GTData).difference({"|","/"}))
                if nHaplotypes>1:
                    allHetIndexes.append(index)
            index+=1

        SNPVCFFile.close()

        random.shuffle(allHetIndexes)
        allAllowedHetIndexes=allHetIndexes[0:nHetAllowed]
        allAllowedHetIndexes.sort()

        i=0
        for line in allLines:
            if i not in allHetIndexes or i in allAllowedHetIndexes:
                hetSNPVCFText+=line
            i+=1

        hetSNPVCFFile=open(subsampledVCFFilePath,"w")
        hetSNPVCFFile.write(hetSNPVCFText)
        hetSNPVCFFile.close()



if __name__=="__main__":
    # Load parameters
    benchmarkParameterFilePath = sys.argv[1]
    runningMode = sys.argv[2] #either "groundTruth" or "virtualHybrids
    benchmarkParameters = loadBenchmarkParameterFile(benchmarkParameterFilePath)

    # Handle paths
    mainPath=benchmarkParameters["mainPath"]
    logPath=mainPath+"log/"
    if runningMode=="groundTruth":
        groundTruthPath = mainPath + "groundTruth/"
        mappedGT=groundTruthPath+"Mapped/"
        allPaths = [groundTruthPath, mappedGT]
    elif runningMode=="virtualHybrids":
        hybridPath = mainPath + "virtualPolyploids/"
        allPaths = [hybridPath]

    if "shortReads" in benchmarkParameters.keys():
        shortReadPath=mainPath+"rawData/shortReads/"
        if runningMode=="groundTruth":
            mappedGTShortReads=mappedGT+"shortReads/"
            variantCalledGT=groundTruthPath+"Variants/"
            variantCalledGTShortReads=variantCalledGT+"shortReads/"
            allPaths+=[shortReadPath,mappedGTShortReads,variantCalledGT,variantCalledGTShortReads]
        elif runningMode=="virtualHybrids":
            shortReadPath=hybridPath+"rawReads/shortReads/"
            mappedShortReads=hybridPath+"Mapped/shortReads/"
            variantCalledShortReads=hybridPath+"Variants/shortReads/"
            allPaths+=[shortReadPath,mappedShortReads,variantCalledShortReads]

    if "longReads" in benchmarkParameters.keys():
        longReadPath=mainPath+"rawData/longReads/"
        if runningMode=="groundTruth":
            mappedGTLongReads=mappedGT+"longReads/"
            allPaths+=[longReadPath,mappedGTLongReads]
        elif runningMode=="virtualHybrids":
            longReadPath=hybridPath+"rawReads/longReads/"
            mappedLongReads=hybridPath+"Mapped/longReads/"
            allPaths+=[longReadPath,mappedLongReads]

    for path in allPaths:
        os.makedirs(path, exist_ok=True)

    if runningMode=="groundTruth":
        fullLogPath=logPath+"groundTruthReadProcessingLog.txt"
    elif runningMode=="virtualHybrids":
        fullLogPath=logPath+"virtualHybridsReadProcessingLog.txt"

    logText=""

    #Mapping parameters
    speciesName=benchmarkParameters["reference"]["speciesName"]
    referencePath=mainPath+"reference/"+speciesName+"/"
    referenceFilePath=referencePath+speciesName+".fasta"
    splitReadFlag="260" #This flag allows split reads
    threads=benchmarkParameters["threads"]

    if runningMode=="groundTruth" and "longReads" in benchmarkParameters.keys():
        longReadMethod=benchmarkParameters["longReadMethod"]
        longReadStrains = list(benchmarkParameters["longReads"].keys())
        for strainName in longReadStrains:
            longReadName = benchmarkParameters["longReads"][strainName]
            longReadFile = longReadPath + longReadName + ".fastq.gz"
            outputLog, systemMessage = longReadMapping(strainName, longReadFile, referenceFilePath, mappedGTLongReads, splitReadFlag, longReadMethod, threads)
            print(systemMessage)
            updateLog(fullLogPath, outputLog)

    if runningMode=="groundTruth" and "shortReads" in benchmarkParameters.keys():
        shortReadStrains=list(benchmarkParameters["shortReads"].keys())
        #Map short reads
        for strainName in shortReadStrains:
            shortReadName = benchmarkParameters["shortReads"][strainName]
            shortReadFile_R1 = shortReadPath + shortReadName + "_1.fastq.gz"
            shortReadFile_R2 = shortReadPath + shortReadName + "_2.fastq.gz"
            outputLog,systemMessage=shortReadMapping(strainName,shortReadFile_R1,shortReadFile_R2,referenceFilePath,mappedGTShortReads,threads)
            print(systemMessage)
            updateLog(fullLogPath,outputLog)
        #Variant calling short reads
        estimatedPloidy = "2"  #Required for GATK's variant calling and ground truth strains are in theory either haploids or homozygous diploids
        strains = ["ACA", "BMB", "CCN", "CRL"]
        pool = Pool()

        results = []
        for strainName in shortReadStrains:
            result = pool.apply_async(runVariantCalling, [strainName,mappedGTShortReads,variantCalledGTShortReads,referenceFilePath,estimatedPloidy])
            results.append(result)

        for result in results:
            result.get()

        print("Variant calling over")

    if runningMode == "groundTruth":
        print("All ground truth reads processed")

    if runningMode == "virtualHybrids" and "longReads" in benchmarkParameters.keys():
        for strainList in benchmarkParameters["strainLists"]:
            for coverageLevel in benchmarkParameters["coverages"]:
                longReadName="_".join(strainList)+"_"+coverageLevel+"X"
                longReadFile=longReadPath+longReadName+".fastq.gz"
                outputLog,systemMessage=longReadMapping(longReadName,longReadFile,referenceFilePath,mappedLongReads,splitReadFlag,"ont",threads)
                print(systemMessage)
                updateLog(fullLogPath,outputLog)

    #Map short reads
    if runningMode == "virtualHybrids" and "shortReads" in benchmarkParameters.keys():
        for strainList in benchmarkParameters["strainLists"]:
            for coverageLevel in benchmarkParameters["coverages"]:
                shortReadName="_".join(strainList)+"_"+coverageLevel+"X"
                shortReadFile_R1 = shortReadPath + shortReadName + "_1.fastq.gz"
                shortReadFile_R2 = shortReadPath + shortReadName + "_2.fastq.gz"
                outputLog,systemMessage=shortReadMapping(shortReadName,shortReadFile_R1,shortReadFile_R2,referenceFilePath,mappedShortReads,threads)
                print(systemMessage)
                updateLog(fullLogPath,outputLog)

        pool = Pool()

        results=[]

        for strainList in benchmarkParameters["strainLists"]:
            strainName="_".join(strainList)
            estimatedPloidy=str(len(strainList))
            result=pool.apply_async(runVariantCalling,[strainName,mappedShortReads,variantCalledShortReads,referenceFilePath,estimatedPloidy])
            results.append(result)

        for result in results:
            result.get()

        print("Variant calling over")

        for strainList in benchmarkParameters["strainLists"]:
            for coverageLevel in benchmarkParameters["coverages"]:
                strainName = "_".join(strainList)+"_"+coverageLevel+"X"
                estimatedPloidy = str(len(strainList))
                for subsamplingSize in benchmarkParameters["heterozygosityRates"]:
                    VCFFile=os.path.join(variantCalledShortReads,strainName+"_"+estimatedPloidy+"n.vcf")
                    subsampleVCF(VCFFile,subsamplingSize,True,referenceFilePath,benchmarkParameters["genomeSize"])
                    subsampleVCF(VCFFile,subsamplingSize,False,referenceFilePath,benchmarkParameters["genomeSize"])

        print("Subsampling over")

    if runningMode == "virtualHybrids":
        print("All virtual hybrid reads processed")
