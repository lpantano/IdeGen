import subprocess
import os

def createIntrons(inv,data,log):
	gen = ''
	out = open(data['out'] + "/gtf/" + inv + "_introns.bed",'w')
	for line in open(data['out'] + "/gtf/" + inv + "_filter.gtf",'r'):
		cols = line.strip().split("\t")
		ids = cols[8].split(";")
		if cols[2]=="exon":
			tr = ids[1].split(" ")[2].replace('"','')
			if gen == tr:
				idx += 1
				intron = tr + "_i" + str(idx)
				end = cols[3]
				out.write("%s\t%s\t%s\t%s\n" % (cols[0],start,end,intron))
				start = cols[4]
			else:
				gen = tr
				start = cols[4]
				idx = 0
	out.close()

def createIntronsBam(inv,data,log):
	cmd = ' '.join(['samtools view -h ',data['out'] + "/bams/" + inv + ".sort.bam | awk '$6~/N/ || $0~/@/' | samtools view -Sb - |bedtools intersect -abam - -b",data['out'] + "/gtf/" + inv + "_filter.gtf","| bedtools bamtobed -split -i -  "])
	p=subprocess.call(cmd, shell=True,stderr=file(data['tempdir']+"/log.bedtools",'w'),stdout=file(data['out'] + "/gtf/" + inv + ".bed",'w'))
	gen = ''
	out = open(data['out'] + "/gtf/" + inv + "_bam_introns.bed",'w')
	for line in open(data['out'] + "/gtf/" + inv + ".bed",'r'):
		cols = line.strip().split("\t")
		tr = cols[3]
		if gen == tr:
			idx += 1
			intron = tr + "_i" + str(idx)
			end = cols[1]
			out.write("%s\t%s\t%s\t%s\n" % (cols[0],start,end,intron))
			start = cols[2]
		else:
			gen = tr
			start = cols[2]
			idx = 0
	out.close()

def countJuntions(inv,data,log):
	cmd = ' '.join(['bedtools intersect -r -f 0.95 -wo  -a ',data['out'] + "/gtf/" + inv + "_introns.bed",' -b ',data['out'] + "/gtf/" + inv + "_bam_introns.bed | awk '($2-5)<$6 && ($2+5)>$6 && ($3-5)<$7 && ($3+5)>$7'| cut -f 4 | sort | uniq -c "])
	p=subprocess.call(cmd, shell=True,stderr=file(data['tempdir']+"/log.counts",'w'),stdout=file(data['out'] + "/gtf/" + inv + ".junctions",'w'))
	return(p)

def convertGenes(inv,data,log):
	cmd = ' '.join(["$INVFUSION/genesConversion.sh",inv,data['refdir']+"/mask.bed",data['refdir'],data['annotation'],data['tempdir'],data['out']+"/gtf"])
	p=subprocess.call(cmd, shell=True,stderr=file(data['tempdir']+"/log.convert",'aw'),stdout=file(data['tempdir']+"/log.convert",'aw'))
	return(p)

def overlapKnownGenes(inv,data,log):
	afile = data['out']+"/gtf/"+inv+".knowngenes.gtf"
	bfile = data['out']+"/gtf/"+inv+"_filter.gtf"
	cmd = ' '.join(['bedtools coverage -hist -a ',afile,'-b',bfile])
	p=subprocess.call(cmd, shell=True,stderr=file(data['tempdir']+"/log.overlapKnownGenes",'aw'),stdout=file(data['out']+"/gtf/"+inv+".overlap",'w'))
	return(p)

def knownExon(inv,data,log):
	p = convertGenes(inv,data,log)
	p = overlapKnownGenes(inv,data,log)
	return(p)

def complexity(inv,data,log):
	afile = data['out']+"/gtf/"+inv+".knowngenes.gtf"
	cmd = ' '.join(['bedtools intersect -wa -c -a ',afile,'-b',afile])
	p=subprocess.call(cmd, shell=True,stderr=file(data['tempdir']+"/log.overlapKnownGenes",'aw'),stdout=file(data['out']+"/gtf/"+inv+".complexity",'w'))
	return(p)
	
def intronBP(inv,data,log):
	bfile = data['out']+"/gtf/"+inv+"_introns.bed"
	afile = data['refdir']+"/"+inv+".bed"
	cmd = ' '.join(['bedtools coverage -hist -a ',afile,'-b',bfile,"| grep -v all  | awk '$8<0.9'"])
	p=subprocess.call(cmd, shell=True,stderr=file(data['tempdir']+"/log.overlapKnownGenes",'aw'),stdout=file(data['out']+"/gtf/"+inv+".intronBP",'w'))
	return(p)


def getFeatures(inv,data,log):
	p = createIntrons(inv,data,log)
	p = createIntronsBam(inv,data,log)
	p = countJuntions(inv,data,log)
	p = knownExon(inv,data,log)
	p = complexity(inv,data,log)
	p = intronBP(inv,data,log)
	return(0)