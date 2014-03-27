import subprocess
import os

def createIndex(data,log):
	log.info("index with STAR")
	newpath = data['refdir']+"/genome"
	listFasta = ""

	if not os.path.exists(newpath):
		os.makedirs(newpath)

	for item in os.listdir(data['refdir']):
		if item.find(".fa")>0: listFasta = listFasta + " " + data['refdir'] + "/" + item 

	cmd=" ".join(['STAR --runMode genomeGenerate --genomeDir ',newpath,' --genomeFastaFiles',listFasta,' --runThreadN ',str(data['threads'])])
	log.debug(cmd)
	p=subprocess.call(cmd, shell=True,stderr=subprocess.STDOUT,stdout=file(data['tempdir']+"/log.index.genome",'w'))
	log.info(p)
	return(2)

def mapGenome(data,log):
	log.info("mapping with STAR")
	genomepath = data['refdir']+"/genome"
	newpath = data['out']+"/bams"
	if not os.path.exists(newpath): os.makedirs(newpath)
	for l in open(data['fasta'],'r'):
		cols = l.strip().split("\t")
		output = newpath + "/" + cols[0] + ".bam"
		p = 0
		if not os.path.exists(output):
			log.info(cols[0])
			cmd = ' '.join(["STAR  --genomeDir",genomepath,"genome --readFilesIn",cols[1],cols[2],"--readFilesCommand zcat --runThreadN 5 --outStd SAM --outSAMmode Full --outSAMattributes All --outFilterType BySJout --outSJfilterReads Unique   --outFilterMultimapNmax 1 --outSAMstrandField intronMotif  | samtools view -bS - >",output])
			log.debug(cmd)
			p=subprocess.call(cmd, shell=True,stderr=file(data['tempdir']+"/log.map.genome",'w'))
			log.info(p)
	return(p)


def align(data,log):
	e = createIndex(data,log)
	e = mapGenome(data,log)
	return(e)