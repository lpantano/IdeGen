import subprocess
import os
from libs.utils import run


def mergeInvcmd(listindv,inversion,data):
	newpath = data['out']+"/bams"
	outfile = newpath + "/" + inversion + ".bam"
	infile = data['tempdir'] + "/" + "inversion.list"
	inf = open(infile,'w')
	for f in listindv:
		inf.write(newpath+"/"+f+".grep.sort.bam\n")
	inf.close()
	cmd = ' '.join(['bamtools filter -region',inversion,' -list',infile,'-out',outfile])
	p=run('mergeInvcmd',cmd, data['tempdir']+"/log.merge",data['tempdir']+"/log.merge")
	return(p)

def sortInvcmd(inversion,data):
	newpath = data['out']+"/bams"
	infile = newpath + "/" + inversion + ".bam"
	outfile = newpath + "/" + inversion + ".sort"
	cmd = ' '.join(['samtools sort ',infile,outfile])
	p=run('sortInvcmd',cmd, data['tempdir']+"/log.sort")
	return(p)

def grepBamcmd(infile,outfile,data):
	cmd = ' '.join(['bamtools filter -script $INVFUSION/libs/filterbam -in',infile,'-out',outfile])
	p=run('grepBamcmd',cmd, data['tempdir']+"/log.grep")
	return(p)

def sortBamcmd(infile,outfile,data):
	cmd = ' '.join(['samtools sort ',infile,outfile])
	p=subprocess.call(cmd, shell=True,stderr=file(data['tempdir']+"/log.sort",'w'))
	return(p)

def indexBamcmd(infile,data):
	cmd = ' '.join(['samtools index ',infile])
	p=run('indexBamcmd',cmd, data['tempdir']+"/log.index")
	return(p)


def sortBam(data,log):
	log.info("Sort bams")
	newpath = data['out']+"/bams"
	if not os.path.exists(newpath): os.makedirs(newpath)
	for l in open(data['fasta'],'r'):
		cols = l.strip().split("\t")
		log.info(cols[0])
		infile = newpath + "/" + cols[0] + ".bam"
		outfile = newpath + "/" + cols[0] + ".sort.bam"
		if not os.path.exists(outfile):
			cmd = ' '.join(['bamtools sort -in',infile,'-out',outfile])
			p=run('sortBam with'+infile,cmd, data['tempdir']+"/log.sort")
			log.info(p)
			if p == 0:
				#p = os.remove(infile)
				log.info("remove file")
	log.info(p)
	return(p)

def prepareBam(data,log):
	log.info("Grep bams")
	newpath = data['out']+"/bams"
	p = 0
	for l in open(data['fasta'],'r'):
		cols = l.strip().split("\t")
		infile = newpath + "/" + cols[0] + ".bam"
		outfile = newpath + "/" + cols[0] + ".grep.bam"
		if not os.path.exists(outfile):
			log.info(cols[0])
			p = grepBamcmd(infile,outfile,data)
		outfile = newpath + "/" + cols[0] + ".grep.sort"
		infile = newpath + "/" + cols[0] + ".grep.bam"
		if not os.path.exists(outfile+".bam"):
			p = sortBamcmd(infile,outfile,data)
		outfile = newpath + "/" + cols[0] + ".grep.sort.bam"
		if not os.path.exists(outfile+".bai"):
			p = indexBamcmd(outfile,data)
		log.info(p)
	return(p)

def getGroup(ma):
	f = open(ma,'r')
	h = f.readline().strip().split("\t")
	hl = {}
	hname = []
	lenl = len(h)
	for hc in h:
		hl[hc] = []
		hname.append(hc)	
	for l in f:
		h = l.strip().split("\t")
		i = h[0]
		#print(i)
		for c in range(1,lenl):
			#print(h[c])
			if h[c] == "INV" :
				hl[hname[c]].append(i) 		
	return(hl)

def CLASScmd(inv,listindv,data,log):
	f = open(data['list'],'r')
	p = 0
	#merge individuals
	p = mergeInvcmd(listindv,inv,data)
	#sort
	p = sortInvcmd(inv,data)
	#run class
	log.info("Running CLASS")
	newpath = data['out']
	infile = newpath + "/bams/" + inv + ".sort" 
	output = newpath + "/gtf/" + inv + ".gtf"
	cmd = ' '.join(['$INVFUSION/class.sh',infile,output])
	log.debug(cmd)
	p = run('CLASScmd',cmd, data['tempdir']+"/log.class")
	#log.info(p)
	return(p)


def launch(data,log):
	f = open(data['list'],'r')
	p = 0
	il = getGroup(data['genotype'])
	newpath = data['out']+"/gtf"
	if not os.path.exists(newpath): os.makedirs(newpath)
	#log.info(il)
	for line in f:
		cols = line.strip().split("\t")
		log.info(cols[0])
		log.info(len(il[cols[0]]))
		if not os.path.exists(newpath + "/" + cols[0] + ".gtf"):
			#just if there are two groups
			p = CLASScmd(cols[0],il[cols[0]],data,log)
			log.info(p)
		#return(0)
	return(p)

def createTranscripts(data,log):
	log.info("de novo transcript assembly")
	p = 0
	p = prepareBam(data,log)
	if p == 0:
		p = launch(data,log)
	return(p)
