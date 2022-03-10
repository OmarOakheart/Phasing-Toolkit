import sys
import os
import random

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

def loadVCFInformation(vcfFilePath):
    VCFInformation={}
    vcfFile=open(vcfFilePath,"r")
    for line in vcfFile:
        line=line.strip("\n")
        if "#" not in line:
            line=line.split("\t")
            chromosome=line[0]
            position=line[1]
            SNPs=[line[3]]+line[4].split(",")
            VCFInformation.setdefault(chromosome,{})
            VCFInformation[chromosome][position]=SNPs
    vcfFile.close()
    return VCFInformation

def floppResultTranslator(floppResultFilePath,outputFilePath,VCFInformation):
    translatedFloppResults=""
    floppFile=open(floppResultFilePath,"r")
    chromosomeName=None
    haplotigIndex=0
    haplotigNames=set()
    for line in floppFile:
        line=line.strip("\n")
        if line=="*****":
            haplotigIndex=max(haplotigNames)
        elif line[0:2]=="**" and line[-2:]=="**":
            chromosomeName=line[2:-2]
        else:
            line=line.split("\t")
            SNPPosition=line[0].split(":")[1]
            SNPCodes=line[1:1+(len(line)-1)//2] #There will always be position+nSNPs+nSNPs in the results
            localIndex=0
            for SNPCode in SNPCodes:
                localIndex += 1
                currentHaplotigName = "flopp_" + str(haplotigIndex + localIndex)
                haplotigNames.add(haplotigIndex + localIndex)
                if int(SNPCode)==-1:
                    SNP="*" #-1 means no read coverage, which means indel.
                else:
                    SNP=VCFInformation[chromosomeName][SNPPosition][int(SNPCode)]
                translatedFloppResults+="\t".join([currentHaplotigName,chromosomeName,SNPPosition,SNP])+"\n"
    floppFile.close()
    outputFile=open(outputFilePath,"w")
    outputFile.write(translatedFloppResults)
    outputFile.close()

def whatsHapPolyphaseResultTranslator(whatsHapPolyphaseFilePath,outputFilePath):
    translatedWhatsHapPolyphaseResults=""
    whatsHapPolyphaseFile=open(whatsHapPolyphaseFilePath,"r")
    for line in whatsHapPolyphaseFile:
        line=line.strip("\n")
        if line[0]=="#":
            pass
        else:
            line=line.split("\t")
            if "PS" in line[8] and "HS" in line[8]:
                chromosomeName=line[0]
                SNPPosition=line[1]
                SNPs=[line[3]]+line[4].split(",")
                infoBlock=line[9]
                GTBlock=infoBlock.split(":")[0]
                phaseBlocks=infoBlock.split(":")[-1].split(",")
                localPhase=0
                for SNPCode in GTBlock.split("|"):
                    localPhase+=1
                    currentHaplotigName=chromosomeName+"_"+phaseBlocks[localPhase-1]+"_"+str(localPhase)
                    SNP=SNPs[int(SNPCode)]
                    translatedWhatsHapPolyphaseResults+="\t".join([currentHaplotigName,chromosomeName,SNPPosition,SNP])+"\n"
    whatsHapPolyphaseFile.close()
    outputFile=open(outputFilePath,"w")
    outputFile.write(translatedWhatsHapPolyphaseResults)
    outputFile.close()

#Now a function to generate the "true" phases using the original VCFs
def getTruePhases(originalVCFs):
    referenceDict={}
    truePhaseDict={}
    heterozygosityChecker={}
    sampleNames=[]
    for originalVCFFilePath, sampleName in originalVCFs:
        sampleNames.append(sampleName)
        originalVCFFile=open(originalVCFFilePath,"r")
        for line in originalVCFFile:
            line=line.strip("\n")
            if line[0]=="#":
                pass
            else:
                line=line.split("\t")
                chromosome=line[0]
                SNPPosition=line[1]
                refSNP=line[3]
                SNPs=line[4].split(",")
                for SNP in SNPs:
                    referenceDict.setdefault(chromosome,{})
                    referenceDict[chromosome][SNPPosition]=chromosome+":"+SNPPosition+"="+refSNP
                    truePhaseDict.setdefault(chromosome,{})
                    truePhaseDict[chromosome].setdefault(SNPPosition, {})
                    truePhaseDict[chromosome][SNPPosition].setdefault(sampleName,set())
                    fullSNP=chromosome+":"+SNPPosition+"="+SNP
                    truePhaseDict[chromosome][SNPPosition][sampleName].add(fullSNP)
                    heterozygosityChecker.setdefault(chromosome, {})
                    heterozygosityChecker[chromosome].setdefault(SNPPosition, set())
                    heterozygosityChecker[chromosome][SNPPosition].add(SNP)
        originalVCFFile.close()
    #Must include the ref allele for positions present in one of the VCFs but not others (Keep track of ref)
    for chromosome, SNPPositions in referenceDict.items():
        for SNPPosition, refHaplotype in SNPPositions.items():
            for sampleName in sampleNames:
                if sampleName not in truePhaseDict[chromosome][SNPPosition].keys():
                    truePhaseDict[chromosome][SNPPosition][sampleName]={refHaplotype}
                    heterozygosityChecker[chromosome][SNPPosition].add(refHaplotype)
    #Only want to include heterozygous SNPs
    heterozygousPhaseDict={}
    for chromosome, chromosomeData in truePhaseDict.items():
        heterozygousPhaseDict.setdefault(chromosome,{})
        for SNPPosition, haplotypeData in chromosomeData.items():
            for sampleName, sampleSNPs in haplotypeData.items():
                heterozygousPhaseDict[chromosome].setdefault(sampleName,set())
                for fullSNP in sampleSNPs:
                    if len(heterozygosityChecker[chromosome][SNPPosition])>1:
                        heterozygousPhaseDict[chromosome][sampleName].add(fullSNP)
    return heterozygousPhaseDict

#Now a function to load translated results
def loadPhasingResults(phasingResultFilePath):
    phasingResults={}
    phasingResultFile=open(phasingResultFilePath,"r")
    for line in phasingResultFile:
        line=line.strip("\n").split("\t")
        haplotigName=line[0]
        chromosome=line[1]
        position=line[2]
        SNP=line[3]
        fullSNP=chromosome+":"+position+"="+SNP
        phasingResults.setdefault(chromosome, {})
        phasingResults[chromosome].setdefault(haplotigName,set())
        phasingResults[chromosome][haplotigName].add(fullSNP)
    phasingResultFile.close()
    return phasingResults

#Now a function to compare & calculate, when given a standardized file
def calculateAccuracyMetrics(groundTruthInfo,predictedInfo,testInfo,accuracyMetricFile):
    accuracyData={"bestScore":{"Total":0},"totalTruePositives":{},"totalFalsePositives":{},"totalFalseNegatives":{},"results":{}}
    #Initializing
    haplotypeNames=set()
    nContigs=0
    for chromosome, haplotypes in groundTruthInfo.items():
        for haplotypeName, SNPs in haplotypes.items():
            nContigs+=1
            haplotypeNames.add(haplotypeName)
            accuracyData["bestScore"].setdefault(haplotypeName,0)
            accuracyData["bestScore"][haplotypeName]+=len(SNPs)
            accuracyData["bestScore"]["Total"]+=len(SNPs)
            accuracyData["totalTruePositives"].setdefault(haplotypeName, 0)
            accuracyData["totalFalsePositives"].setdefault(haplotypeName,0)
            accuracyData["totalFalseNegatives"].setdefault(haplotypeName,set())
            SNPPos=set([fullSNP.split("=")[0] for fullSNP in SNPs])
            accuracyData["totalFalseNegatives"][haplotypeName]=accuracyData["totalFalseNegatives"][haplotypeName].union(SNPPos)

    nHaplotigs=0
    for chromosome, haplotypes in predictedInfo.items():
        for predictionName, SNPs in haplotypes.items():
            nHaplotigs+=1
            closestHaplotypeName,TP,FP,FNPos=getBestMatch(SNPs,groundTruthInfo)
            if closestHaplotypeName=="None":
                closestHaplotypeName=random.choice(list(haplotypeNames))
            accuracyData["totalTruePositives"][closestHaplotypeName]+=TP
            accuracyData["totalFalsePositives"][closestHaplotypeName]+=FP
            accuracyData["totalFalseNegatives"][closestHaplotypeName].difference_update(FNPos)

    P90Value=getPn(predictedInfo,groundTruthInfo,90)
    accuracyData["results"]["P90"]=P90Value
    accuracyData["results"]["nHaplotigs"]=nHaplotigs
    accuracyData["results"]["nContigs"]=nContigs
    fullTP=0
    fullFP=0
    fullFN=0
    for haplotypeName in haplotypeNames:
        accuracyData["results"].setdefault(haplotypeName,{"TP":0,"FP":0,"FN":0})
        TP=accuracyData["totalTruePositives"][haplotypeName]
        FP=accuracyData["totalFalsePositives"][haplotypeName]
        FN=len(accuracyData["totalFalseNegatives"][haplotypeName])
        total=TP+FP+FN
        fullTP+=TP
        fullFP+=FP
        fullFN+=FN
        accuracyData["totalFalseNegatives"][haplotypeName]=len(accuracyData["totalFalseNegatives"][haplotypeName])
        accuracyData["results"][haplotypeName]["TP"]=round(TP/total*100,2)
        accuracyData["results"][haplotypeName]["FP"]=round(FP/total*100,2)
        accuracyData["results"][haplotypeName]["FN"]=round(FN/total*100,2)
    fullTotal=fullTP+fullFP+fullFN
    accuracyData["results"]["all"]={}
    accuracyData["results"]["all"]["TP"]=round(fullTP/fullTotal*100,2)
    accuracyData["results"]["all"]["FP"]=round(fullFP/fullTotal*100,2)
    accuracyData["results"]["all"]["FN"]=round(fullFN/fullTotal*100,2)
    resultTable=testInfo+[accuracyData["results"]["all"]["TP"],accuracyData["results"]["all"]["FP"],accuracyData["results"]["all"]["FN"],accuracyData["results"]["P90"],accuracyData["results"]["nHaplotigs"],accuracyData["results"]["nContigs"]]
    accuracyLine="\t".join([str(element) for element in resultTable])
    pFile = open(accuracyMetricFile, "a")
    pFile.write(accuracyLine+"\n")
    pFile.close()
    print(accuracyLine)

def getPn(predictedInfo,groundTruthInfo,n):
    totalSNPs={"Total":0}
    SNPsToPredict={}
    for chromosome, haplotypes in groundTruthInfo.items():
        for haplotypeName, SNPs in haplotypes.items():
            SNPPos=set([fullSNP.split("=")[0] for fullSNP in SNPs])
            totalSNPs.setdefault(haplotypeName,0)
            totalSNPs[haplotypeName]+=len(SNPPos)
            totalSNPs["Total"]+=len(SNPPos)
            SNPsToPredict.setdefault(haplotypeName,set())
            SNPsToPredict[haplotypeName]=SNPsToPredict[haplotypeName].union(SNPPos)

    contributors=[]
    for chromosome, haplotypes in predictedInfo.items():
        for predictionName, SNPs in haplotypes.items():
            closestHaplotypeName,TP,FP,predictedSNPs=getBestMatch(SNPs,groundTruthInfo)
            contributors.append((len(predictedSNPs),predictionName,closestHaplotypeName,predictedSNPs))

    contributors.sort(reverse=True)

    Pn=0
    for contributor in contributors:
        Pn+=1
        haplotypeName=contributor[2]
        predictedSNPs=contributor[3]
        if haplotypeName=="None":
            pass
        else:
            SNPsToPredict[haplotypeName].difference_update(predictedSNPs)
        remainingScore=0
        for haplotype, haplotypeSNPPositions in SNPsToPredict.items():
            remainingScore+=len(haplotypeSNPPositions)
        if (1-remainingScore/totalSNPs["Total"])*100>n:
            return Pn
    return -1
    #This gets you Pn, for example the P90, which is the minimum number of contigs needed to phase 90% of SNPs once.

def getBestMatch(predictedSNPs,groundTruthInfo):
    bestScore=0
    bestScoreDetail=("None",0,0,set())
    for SNP in predictedSNPs:
        chromosome=SNP.split(":")[0]
        break
    for haplotypeName, trueSNPs in groundTruthInfo[chromosome].items():
        TPScore=len(predictedSNPs&trueSNPs)
        posSNPs=set([fullSNP.split("=")[0] for fullSNP in predictedSNPs])
        posTrueSNPs=set([fullSNP.split("=")[0] for fullSNP in trueSNPs])
        commonPos=posSNPs&posTrueSNPs
        FPScore=0
        FPCandidates=predictedSNPs.difference(trueSNPs)
        for FPCandidate in FPCandidates:
            if FPCandidate.split("=")[0] in commonPos:
                FPScore+=1
        commonPosNumber=len(commonPos)
        noiseScore = len(posSNPs.difference(posTrueSNPs))  # Sometimes there are predictions that aren't in the ground truth, not sure how to handle this
        missingScore=commonPosNumber+noiseScore-(TPScore+FPScore) #Adding noiseScore since it's not a TP, or a FP in the sense we're using it, or missing.
        localScore=TPScore/(TPScore+FPScore+missingScore)
        if localScore>bestScore:
            bestScore=localScore
            bestScoreDetail=(haplotypeName,TPScore,FPScore,commonPos)
    return bestScoreDetail

def subsetTruePhases(groundTruthDict,subsetFile,hybridVCFPrefix):
    VCFInfo=loadVCFInformation(hybridVCFPrefix+subsetFile)
    subsetPositions=set()
    for chromosome, chromosomeData in VCFInfo.items():
        for SNPPosition in chromosomeData.keys():
            subsetPositions.add(chromosome+":"+SNPPosition)
    truePhaseSubset={}
    for chromosome, haplotypeData in groundTruthDict.items():
        truePhaseSubset.setdefault(chromosome,{})
        for haplotypeName, SNPSet in haplotypeData.items():
            truePhaseSubset[chromosome].setdefault(haplotypeName,set())
            for fullSNP in SNPSet:
                if fullSNP.split("=")[0] in subsetPositions:
                    truePhaseSubset[chromosome][haplotypeName].add(fullSNP)
    return truePhaseSubset

if __name__ == "__main__":
    # Load parameters
    benchmarkParameterFilePath = sys.argv[1]
    benchmarkParameters = loadBenchmarkParameterFile(benchmarkParameterFilePath)
    phasingMethods=sys.argv[2:]

    # Handle paths
    mainPath=benchmarkParameters["mainPath"]
    accuracyMetricPath=mainPath+"accuracyMetrics/"
    accuracyMetricFile=accuracyMetricPath+"accuracyMetrics.tsv"
    VCFPrefix=mainPath+"groundTruth/Variants/shortReads/"
    accuracyCalculationPath=mainPath+"accuracyCalculations/"
    WHPPrefix=mainPath+"/phasingPredictions/whatsHapPolyphase/"
    floppPrefix=mainPath+"/phasingPredictions/flopp/"
    nPhasePrefix=mainPath+"/phasingPredictions/nPhase/"
    hybridVCFPrefix=mainPath+"/virtualPolyploids/Variants/shortReads/"

    allPaths=[accuracyMetricPath,accuracyCalculationPath]

    logPath=mainPath+"log/"

    for path in allPaths:
        os.makedirs(path, exist_ok=True)

    #Prepare log

    fullLogPath=logPath+"phasingLog.txt"
    logText=""

    pFile = open(accuracyMetricFile, "a")
    pFile.write("#tool\tploidy\theterozygosity level (%)\tcoverage\tindelStatus\tTrue Positives (%)\tFalse Positives (%)\tmissing (%)\tP90\tnHaplotigs\tnContigs\n")
    pFile.close()

    truePhaseDict={} #Will hold the ground truth of individual virtual strains
    allSubsets={} #Will hold the subset ground truth of individual virtual strains
    for strainList in benchmarkParameters["strainLists"]:
        virtualStrainName="_".join(strainList)
        truePhaseGenotypes=[]
        for strainName in strainList:
            strainVCF=VCFPrefix+strainName+".vcf"
            truePhaseGenotypes.append((strainVCF,strainName))
        truePhaseDict[virtualStrainName]=getTruePhases(truePhaseGenotypes) #The phases we get here are full phases but we generate subsets
        #So we want it to only consider the subsets
        for heterozygosityRate in benchmarkParameters["heterozygosityRates"]:
            for coverageLevel in benchmarkParameters["coverages"]:
                subsetFileName=virtualStrainName+"_"+coverageLevel+"X_"+heterozygosityRate+".SNPs.vcf"
                subsetFileNameIndels=virtualStrainName+"_"+coverageLevel+"X_"+heterozygosityRate+".vcf"
                allSubsets[subsetFileName]=subsetTruePhases(truePhaseDict[virtualStrainName],subsetFileName,hybridVCFPrefix)
                allSubsets[subsetFileNameIndels]=subsetTruePhases(truePhaseDict[virtualStrainName],subsetFileNameIndels,hybridVCFPrefix)

    ################################################################
    #How to pre-process some results & load them in a usable format#
    ################################################################

    if "whatshap-polyphase" in phasingMethods:
        #Translate WH-PP results

        for strainList in benchmarkParameters["strainLists"]:
            virtualStrainName = "_".join(strainList)
            for heterozygosityRate in benchmarkParameters["heterozygosityRates"]:
                for coverage in benchmarkParameters["coverages"]:
                    VCFFileName = virtualStrainName + "_" + coverage + "X_" + heterozygosityRate + ".SNPs.vcf"
                    VCFFileNameIndels = virtualStrainName + "_" + coverage + "X_" + heterozygosityRate + ".vcf"
                    phasingMethod="WHP-PP"
                    ploidy=str(len(strainList))
                    if ".SNPs." in VCFFileName:
                        indelStatus="noIndels"
                    else:
                        indelStatus="Indels"
                    #First no indels
                    testInfo=[phasingMethod,ploidy,heterozygosityRate,coverage,"noIndels"]
                    whatsHapPolyphaseResultTranslator(WHPPrefix+"WHP_"+VCFFileName, accuracyCalculationPath + "Translated_WHP_"+VCFFileName)
                    WHPPhasingResults=loadPhasingResults(accuracyCalculationPath+"Translated_WHP_"+VCFFileName)
                    calculateAccuracyMetrics(allSubsets[VCFFileName],WHPPhasingResults,testInfo,accuracyMetricFile)
                    #Then with indels
                    testInfo=[phasingMethod,ploidy,heterozygosityRate,coverage,"Indels"]
                    whatsHapPolyphaseResultTranslator(WHPPrefix+"WHP_"+VCFFileNameIndels, accuracyCalculationPath + "Translated_WHP_"+VCFFileNameIndels)
                    WHPPhasingResults=loadPhasingResults(accuracyCalculationPath+"Translated_WHP_"+VCFFileNameIndels)
                    calculateAccuracyMetrics(allSubsets[VCFFileNameIndels],WHPPhasingResults,testInfo,accuracyMetricFile)

        print("Done calculating WhatsHap polyphase accuracy")

    if "flopp" in phasingMethods:
        #Translate all flopp results
        for strainList in benchmarkParameters["strainLists"]:
            virtualStrainName = "_".join(strainList)
            for heterozygosityRate in benchmarkParameters["heterozygosityRates"]:
                for coverage in benchmarkParameters["coverages"]:
                    VCFFileName = virtualStrainName + "_" + coverage + "X_" + heterozygosityRate + ".SNPs.vcf"
                    VCFFileNameIndels = virtualStrainName + "_" + coverage + "X_" + heterozygosityRate + ".vcf"
                    phasingMethod="flopp"
                    ploidy=str(len(strainList))
                    #First no indels
                    testInfo=[phasingMethod,ploidy,heterozygosityRate,coverage,"noIndels"]
                    VCFInformation=loadVCFInformation(hybridVCFPrefix+VCFFileName)
                    floppResultTranslator(floppPrefix+"flopp_"+VCFFileName,accuracyCalculationPath+"Translated_flopp_"+VCFFileName,VCFInformation)
                    floppPredictions=loadPhasingResults(accuracyCalculationPath+"Translated_flopp_"+VCFFileName)
                    calculateAccuracyMetrics(allSubsets[VCFFileName], floppPredictions,testInfo,accuracyMetricFile)
                    #Then with indels
                    testInfo=[phasingMethod,ploidy,heterozygosityRate,coverage,"Indels"]
                    VCFInformation=loadVCFInformation(hybridVCFPrefix+VCFFileNameIndels)
                    floppResultTranslator(floppPrefix+"flopp_"+VCFFileNameIndels,accuracyCalculationPath+"Translated_flopp_"+VCFFileNameIndels,VCFInformation)
                    floppPredictions=loadPhasingResults(accuracyCalculationPath+"Translated_flopp_"+VCFFileNameIndels)
                    calculateAccuracyMetrics(allSubsets[VCFFileNameIndels], floppPredictions,testInfo,accuracyMetricFile)
        print("Done calculating flopp accuracy")

    if "nphase" in phasingMethods:
        nPhaseDefaultSuffix="_0.1_0.01_0.05_0_variants.tsv"
        nPhaseDefaultSuffixIndel="_Indels_0.1_0.01_0.05_0_variants.tsv"
        for strainList in benchmarkParameters["strainLists"]:
            virtualStrainName = "_".join(strainList)
            for heterozygosityRate in benchmarkParameters["heterozygosityRates"]:
                for coverage in benchmarkParameters["coverages"]:
                    exactTestName=virtualStrainName + "_" + coverage + "X_" + heterozygosityRate
                    VCFFileName = exactTestName+nPhaseDefaultSuffix
                    VCFFileNameIndels = exactTestName +nPhaseDefaultSuffixIndel
                    phasingMethod="nPhase"
                    ploidy=str(len(strainList))
                    #Without indels
                    testInfo=[phasingMethod,ploidy,heterozygosityRate,coverage,"noIndels"]
                    nPhaseResultPath = nPhasePrefix + exactTestName+"/Phased/" + VCFFileName
                    subsetFileName=exactTestName+".SNPs.vcf"
                    nPhasePhasingResults = loadPhasingResults(nPhaseResultPath)
                    calculateAccuracyMetrics(allSubsets[subsetFileName],nPhasePhasingResults,testInfo,accuracyMetricFile)
                    #With indels
                    testInfo=[phasingMethod,ploidy,heterozygosityRate,coverage,"Indels"]
                    nPhaseIndelResultPath=nPhasePrefix+exactTestName+"_Indels/Phased/"+VCFFileNameIndels
                    subsetFileNameIndels=exactTestName+".vcf"
                    nPhasePhasingResults=loadPhasingResults(nPhaseIndelResultPath)
                    calculateAccuracyMetrics(allSubsets[subsetFileNameIndels],nPhasePhasingResults,testInfo,accuracyMetricFile)

        print("Done calculating nPhase accuracy")

    print("Done calculating accuracy metrics")

