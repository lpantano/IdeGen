import subprocess
import os
from mcmc import estimatemcmc
from toolsGenomic import getFeatures
from candidateClass import candidate
from libs.utils import run

def countscmd(data,gtf,outfile,bamfiles,threads=1,metafeature="-g exon_id -O "):
	cmd = ' '.join(['featureCounts ',metafeature,'-T',str(threads),' -a ',gtf,'-o',outfile,bamfiles])
	p=run('countscmd',cmd, data['tempdir']+"/log.counts")
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
		p=run('normalizationref',cmd,data['tempdir']+"/log.norm",data['tempdir']+"/log.norm")
	return(p)

def normalizationcandidates(inv,data,log):
	infile = data['outcounts'] + "/" + inv + ".tab"
	sizefactors =   data['outcounts'] + "/sizefactors"
	cmd = ' '.join(['Rscript $INVFUSION/R/candidatesnorm.R ',infile,data['genotype'],sizefactors,data['fasta'],inv])
	p = 0
	if os.path.exists(infile):
		if not os.path.exists(infile + ".norm1" ):
			p=run('normalizationcandidates',cmd, data['tempdir']+"/log.cannorm",data['tempdir']+"/log.cannorm")
	return(p)

def mcmc(inv,data,log):
	if os.path.exists(data['outcounts'] + "/" + inv + ".tab.norm1"):
		infile1 = open(data['outcounts'] + "/" + inv + ".tab.norm1",'r')
		infile2 = open(data['outcounts'] + "/" + inv + ".tab.norm2",'r')
		if not os.path.exists(data['outcounts'] + "/" + inv + ".estimator"):
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

def createResults(inv,data,log):
	listres = {}
	if os.path.exists(data['outcounts'] + "/" + inv + ".estimator"):
		for line in open(data['out']+"/gtf/"+inv+"_filter.gtf",'r'):
			cols = line.strip().split("\t")
			if cols[2] == "exon":
				idx=cols[8].split(";")[4].replace('"','').split(" ")[2]
				listres[idx] = candidate(idx)
				listres[idx].tr = cols[8].split(";")[1].replace('"','').split(" ")[2]
		for line in open(data['outcounts'] + "/" +inv+".estimator",'r'):
			cols = line.strip().split("\t")
			idx = cols[0]
			listres[idx].g1=map(float,cols[1:4])
			listres[idx].g2=map(float,cols[4:])
			v = listres[idx].setcf()
		for line in open(data['out']+"/gtf/"+inv+".overlap",'r'):
			cols = line.strip().split("\t")
			if cols[2] == "exon":
				idx=cols[8].split(";")[4].replace('"','').split(" ")[2]
				listres[idx].known=float(cols[12])
		for line in open(data['out']+"/gtf/"+inv+".complexity",'r'):
			cols = line.strip().split("\t")
			if cols[2] == "exon" and line.find("exon_id")>0:			
				idx=cols[8].split(";")[4].replace('"','').split(" ")[2]
				listres[idx].cpx=int(cols[9])-1			
		for line in open(data['out']+"/gtf/"+inv+".intronBP",'r'):
			cols = line.strip().split("\t")
			tr=cols[3].split("_")[0]
			inum=int(cols[3].split("_")[1].replace("i",''))
			idxl = tr+"_"+str(inum)
			idxr = tr+"_"+str(inum+1)
			if listres.has_key(idxl):
				listres[idxl].i2bp = "YES"
			if listres.has_key(idxr):
				listres[idxr].i1bp = "YES"
		for line in open(data['out']+"/gtf/"+inv+".junctions",'r'):
			cols = line.strip().split(" ")
			tr=cols[1].split("_")[0]
			inum=int(cols[1].split("_")[1].replace("i",''))
			idxl = tr+"_"+str(inum)
			idxr = tr+"_"+str(inum+1)
			if listres.has_key(idxl):
				listres[idxl].i2reads = int(cols[0])
			if listres.has_key(idxr):
				listres[idxr].i1reads = int(cols[0])
		indiv = len(open(data['out']+"/bayesian/"+inv+".tab.norm1",'r').readline().strip().split("\t"))-1
		out = open(data['out']+"/svdgResults.tab",'aw')
		outraw = open(data['out']+"/svdg.rda",'aw')
		for idx in listres.keys():
			i = listres[idx]
			sc = 0 
			aveintron1 = 1.0*i.i1reads/(1.*indiv)
			aveintron2 = 1.0*i.i2reads/(1.*indiv)

			outraw.write('\t'.join(map(str,[i.tr,i.number,"exon",i.known]+i.g1+i.g2+["NA"]))+"\n")
			outraw.write('\t'.join(map(str,[i.tr,i.number,"intron",aveintron2]+[0,0,0]+[0,0,0]+[i.i2bp]))+"\n")
			if (i.i1bp=="YES" or i.i2bp=="YES"):
				sc += 1
			if 1.0*i.i1reads/indiv > 3 or .0*i.i2reads/indiv>3:
				sc += 1
			if i.cf >1.8 and i.known<0.1:
				sc += 1
			if  i.cpx <= 3:
				sc += 1
			out.write("%s %s %s %s %s %s %s %s %s %s\n" % (idx,i.tr,i.cf,i.known,i.cpx,i.i1bp,aveintron1,i.i2bp,aveintron2,sc))
		out.close()
	return(0)

def createFigures(data,log):
	newpath = data['out'] + "/img"
	if not os.path.exists(newpath):
		os.makedirs(newpath)
	infile=data['out']+"/svdg.rda"
	cmd = ' '.join(['Rscript $INVFUSION/R/figures.R ',infile,newpath])
	p = 0
	if os.path.exists(infile):
		p=run('createFigures',cmd, data['tempdir']+"/log.fig",data['tempdir']+"/log.fig")
	return(p)


def bayesianEstiamtor(data,log):
	data['outcounts'] = data['out'] + "/bayesian"  
	newpath = data['outcounts']
	if not os.path.exists(newpath): os.makedirs(newpath)
	log.info("featureCounts on reference genes")
	p = countsref(data,log)
	log.info("Calculating size factor")
	p = normalizationref(data,log)
	if os.path.exists(data['out']+"/svdgResults.tab"): os.remove(data['out']+"/svdgResults.tab")
	if os.path.exists(data['out']+"/svdg.rda"): os.remove(data['out']+"/svdg.rda")
	for line in open(data['list'],'r'):
		cols = line.strip().split("\t")
		log.info(cols[0])
		gtf = data['out'] + "/gtf/" + cols[0] + "_filter.gtf"
		statinfo = os.stat(gtf)
		p = 0
		if statinfo.st_size>0:
			log.info("Retrieve genomic features")
			p = getFeatures(cols[0],data,log)
			log.info("Counts canditates")
			p = countscandidates(cols[0],data,log)
			log.info("Normalize candidates")
			p = normalizationcandidates(cols[0],data,log)
			log.info("mcmc")
			p = mcmc(cols[0],data,log)
			log.info("Write results")
			p = createResults(cols[0],data,log)
			#return(0)
			#mcmc
	log.info("Create figures")
	p = createFigures(data,log)
	return(p)