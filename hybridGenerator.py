import sys
import os
import random
import gzip

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

def downsample(coverage,genomeSize,strainName,outputFolder,fastQFiles,readType):
    if readType=="shortReads":
        R1FileName=outputFolder+"/"+strainName+"_1.fastq.gz"
        R2FileName=outputFolder+"/"+strainName+"_2.fastq.gz"
        if os.path.exists(R1FileName)==False and os.path.exists(R2FileName)==False:
            fastQFile1=fastQFiles[0]
            fastQFile2=fastQFiles[1]
            shortReads={"/1":{},"/2":{}}
            i=1
            file1=gzip.open(fastQFile1,"rb")
            for line in file1:
                line=str(line,encoding="utf-8").strip("\n")
                if i%4==1:
                    ID=line.split(" ")[0]
                    if ID[-2:]=="/1" or ID[-2:]=="/2":
                        ID=ID[:-2]
                    RNum="/1"
                    shortReads[RNum][ID]=[]
                else:
                    shortReads[RNum][ID].append(line)
                i+=1
                i=i%4
            file1.close()
            i=1
            file2=gzip.open(fastQFile2,"rb")
            for line in file2:
                line=str(line,encoding="utf-8").strip("\n")
                if i%4==1:
                    ID=line.split(" ")[0]
                    if ID[-2:]=="/1" or ID[-2:]=="/2":
                        ID=ID[:-2]
                    RNum="/2"
                    shortReads[RNum][ID]=[]
                else:
                    shortReads[RNum][ID].append(line)
                i+=1
                i=i%4
            file2.close()

            allReadNames=list(shortReads["/1"].keys())
            random.shuffle(allReadNames)

            newR1File=""
            newR2File=""

            totalBases=0

            print("Downsampling reads")
            for read in allReadNames:
                if totalBases/genomeSize>coverage:
                    break
                R1=shortReads["/1"][read][0]
                R2=shortReads["/2"][read][0]
                newLine1=read+"\n"+"\n".join(shortReads["/1"][read])+"\n"
                newR1File+=newLine1
                newLine2=read+"\n"+"\n".join(shortReads["/2"][read])+"\n"
                newR2File+=newLine2
                totalBases+=len(R1)+len(R2)

            print("Starting to encode")
            newR1File=str.encode(newR1File)
            print("Done encoding R1 reads")
            newR2File=str.encode(newR2File)
            print("Done encoding R2 reads")
            output1=gzip.open(R1FileName,"wb")
            output1.write(newR1File)
            output1.close()
            output2=gzip.open(R2FileName,"wb")
            output2.write(newR2File)
            output2.close()
        else:
            print()

    elif readType=="longReads":
        longReadFile=outputFolder+"/"+strainName+".fastq.gz"
        if os.path.exists(longReadFile) == False:
            fastQFile1=fastQFiles[0]
            longReads={}
            i=1
            file=gzip.open(fastQFile1,"rb")
            for line in file:
                line=str(line,encoding="utf-8").strip("\n")
                if i%4==1:
                        ID=line
                        longReads[ID]=[]
                else:
                        longReads[ID].append(line)
                i+=1
                i=i%4
            file.close()

            allReadNames=list(longReads.keys())
            random.shuffle(allReadNames)

            newFile=""
            totalBases=0

            print("Starting to downsample")
            for read in allReadNames:
                if totalBases/genomeSize>coverage:
                        break
                if read in longReads:
                    try:
                        lineLength=len(longReads[read][0])
                        newLine=read+"\n"+"\n".join(longReads[read])+"\n"
                        newFile+=newLine
                        totalBases+=lineLength
                    except:
                        pass
            print("Starting to encode")
            newFile=str.encode(newFile)
            print("Done encoding")
            output=gzip.open(longReadFile,"wb")
            output.write(newFile)
            output.close()
        else:
            print(longReadFile,"already exists.")

def mergeReads(outputFileName,filesToMerge):
    #Load all the files
    if os.path.exists(outputFileName)==False:
        allReads={}
        for mergeFile in filesToMerge:
            i=1
            openFile=gzip.open(mergeFile,"rb")
            for line in openFile:
                line=str(line,encoding="utf-8").strip("\n")
                if i%4==1:
                    ID=line.split(" ")[0]
                    allReads[ID]=[]
                else:
                    allReads[ID].append(line)
                i+=1
                i=i%4
            openFile.close()
        #Randomly merge them back together
        allReadNames = list(allReads.keys())
        allReadNames.sort()

        mergedFileText = ""

        for read in allReadNames:
            newLine = read + "\n" + "\n".join(allReads[read]) + "\n"
            mergedFileText += newLine

        print("Starting to encode")
        mergedFile = str.encode(mergedFileText)
        print("Done encoding reads")
        output = gzip.open(outputFileName, "wb")
        output.write(mergedFile)
        output.close()
    else:
        print(outputFileName,"already exists.")

def hybridGenerator(strains,strainDictionary,coverage):
    hybridName="_".join(strains)+"_"+str(coverage)+"X"
    genomeSize=int(strainDictionary["genomeSize"])
    shortReadsToMerge_R1=[]
    shortReadsToMerge_R2=[]
    longReadsToMerge=[]
    for strainName in strains:
        downsampledName=strainName+"_"+str(coverage)+"X"
        longReadFileName=strainDictionary["longReads"][strainName]
        shortReadFileName=strainDictionary["shortReads"][strainName]
        shortReadFilePath_R1=shortReadPath+shortReadFileName+"_1.fastq.gz"
        shortReadFilePath_R2=shortReadPath+shortReadFileName+"_2.fastq.gz"
        longReadFilePath=longReadPath+longReadFileName+".fastq.gz"
        #Create short read file at coverage
        downsample(coverage, genomeSize, downsampledName, hybridShortReads, [shortReadFilePath_R1,shortReadFilePath_R2], "shortReads")
        shortReadsToMerge_R1.append(hybridShortReads+downsampledName+"_1.fastq.gz")
        shortReadsToMerge_R2.append(hybridShortReads+downsampledName+"_2.fastq.gz")
        #Create long read at coverage
        downsample(coverage, genomeSize, downsampledName, hybridLongReads, [longReadFilePath], "longReads")
        longReadsToMerge.append(hybridLongReads + downsampledName + ".fastq.gz")
    #Now you want to merge these together.
    hybridShortReadName_R1 = hybridShortReads + hybridName + "_1.fastq.gz"
    hybridShortReadName_R2 = hybridShortReads + hybridName + "_2.fastq.gz"
    hybridLongReadName = hybridLongReads + hybridName + ".fastq.gz"
    mergeReads(hybridShortReadName_R1,shortReadsToMerge_R1)
    print("R1 reads merged")
    mergeReads(hybridShortReadName_R2, shortReadsToMerge_R2)
    print("R2 reads merged")
    mergeReads(hybridLongReadName,longReadsToMerge)
    print("Long reads merged")

if __name__=="__main__":
    # Load parameters
    benchmarkParameterFilePath = sys.argv[1]
    benchmarkParameters = loadBenchmarkParameterFile(benchmarkParameterFilePath)

    # Handle paths
    mainPath=benchmarkParameters["mainPath"]

    shortReadPath=mainPath+"rawData/shortReads/"
    longReadPath=mainPath+"rawData/longReads/"
    hybridPath=mainPath+"virtualPolyploids/"
    hybridRawReads=hybridPath+"rawReads/"
    hybridShortReads=hybridRawReads+"shortReads/"
    hybridLongReads=hybridRawReads+"longReads/"

    logPath=mainPath+"log/"

    allPaths=[hybridPath,hybridRawReads,hybridShortReads,hybridLongReads]

    for path in allPaths:
        os.makedirs(path, exist_ok=True)

    #Prepare log

    fullLogPath=logPath+"hybridGeneratorLog.txt"
    logText=""

    #Prepare virtual polyploids at a coverage of coverageX

    for strains in benchmarkParameters["strainLists"]:
        for coverage in benchmarkParameters["coverages"]:
            hybridGenerator(strains,benchmarkParameters,coverage)

