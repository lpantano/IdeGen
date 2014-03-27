import subprocess
import os
from mcmc import estimatemcmc
from toolsGenomic import getFeatures


def countscmd(data,gtf,outfile,bamfiles,threads=1,metafeature="-g exon_id -O "):
	cmd = ' '.join(['featureCounts -p ',metafeature,'-T',str(threads),' -a ',gtf,'-o',outfile,bamfiles])
	p=subprocess.call(cmd, shell=True,stderr=file(data['tempdir']+"/log.counts",'w'))
	return(p)

def countsref(data,log):
	gtf = data['annotation']
	bamlist = []
	for l in open(data['fasta']):
		cols = l.strip().split("\t")
		bamlist.append(data['out'] + "/bams/" + cols[0] + ".bam")
	bamfiles = ' '.join(bamlist)
	if not os.path.exists(data['outcounts'] + "/reference.tab"):
		p = countscmd(data,gtf,data['outcounts'] + "/reference.tab",bamfiles,data['threads'],"")
	return(0)

def countscandidates(inv,data,log):
	gtf = data['out'] + "/gtf/" + inv + "_filter.gtf"
	statinfo = os.stat(gtf)
	p = 0
	if statinfo.st_size>0:
		bamlist = []
		if not os.path.exists(data['outcounts'] + "/" + inv + ".tab"):
			for l in open(data['fasta']):
				cols = l.strip().split("\t")
				bamlist.append(data['out'] + "/bams/" + cols[0] + ".grep.sort.bam")
			bamfiles = ' '.join(bamlist)
			p = countscmd(data,gtf,data['outcounts'] + "/" + inv + ".tab",bamfiles,data['threads'])
	return(p)


def normalizationref(data,log):
	infile = data['outcounts'] + "/reference.tab"  
	cmd = ' '.join(['Rscript $INVFUSION/R/normalization.R ',infile])
	p = 0
	if not os.path.exists(data['outcounts'] + "/sizefactors" ):
		p=subprocess.call(cmd, shell=True,stderr=file(data['tempdir']+"/log.norm"),stdout=file(data['tempdir']+"/log.norm",'a'))
	return(p)

def normalizationcandidates(inv,data,log):
	infile = data['outcounts'] + "/" + inv + ".tab"
	sizefactors =   data['outcounts'] + "/sizefactors"
	cmd = ' '.join(['Rscript $INVFUSION/R/candidatesnorm.R ',infile,data['genotype'],sizefactors,data['fasta'],inv])
	p = 0
	if os.path.exists(infile):
		if not os.path.exists(infile + ".norm1" ):
			p=subprocess.call(cmd, shell=True,stderr=file(data['tempdir']+"/log.cannorm",'aw'),stdout=file(data['tempdir']+"/log.cannorm",'aw'))
	return(p)

def mcmc(inv,data,log):
	if os.path.exists(data['outcounts'] + "/" + inv + ".tab.norm1"):
		infile1 = open(data['outcounts'] + "/" + inv + ".tab.norm1",'r')
		infile2 = open(data['outcounts'] + "/" + inv + ".tab.norm2",'r')
		out = open(data['outcounts'] + "/" + inv + ".estimator",'w')
		for line1 in infile1:
			line2 = infile2.readline()
			cols1 = line1.strip().split("\t")
			cols2 = line2.strip().split("\t")
			gene = cols1[0]
			values1 = (map(int,cols1[1:]))
			values2 = (map(int,cols2[1:]))
			res = estimatemcmc(gene,inv,values1,values2,data)
			out.write(gene+"\t"+'\t'.join(map(str,res))+"\n")
		out.close()
	return(0)

def createList(inv,data,log):
	return(0)

def bayesianEstiamtor(data,log):
	data['outcounts'] = data['out'] + "/bayesian"  
	newpath = data['outcounts']
	if not os.path.exists(newpath): os.makedirs(newpath)
	log.info("featureCounts on reference genes")
	p = countsref(data,log)
	log.info("Calculating size factor")
	p = normalizationref(data,log)
	for line in open(data['list'],'r'):
		cols = line.strip().split("\t")
		log.info(cols[0])
		gtf = data['out'] + "/gtf/" + cols[0] + "_filter.gtf"
		statinfo = os.stat(gtf)
		p = 0
		if statinfo.st_size>0:
			#p = getFeatures(cols[0],data,log)
			p = countscandidates(cols[0],data,log)
			p = normalizationcandidates(cols[0],data,log)
			p = mcmc(cols[0],data,log)
			return(0)
			#mcmc
	return(0)