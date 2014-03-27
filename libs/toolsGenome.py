import subprocess

def insilico(invname,bedfile,fastaref,tmp,output,size,log):
	log.info(invname)
	cmd=" ".join(['$INVFUSION/makeref.insilico.sh ',invname,bedfile,fastaref,size,tmp,output])
	log.debug(cmd)
	p=subprocess.call(cmd, shell=True,stderr=subprocess.STDOUT,stdout=file(tmp+"/log."+invname,'w'))
	log.info(p)
	return(p)

def real(invname,bedfile,genbank,fastaref,tmp,output,size,log):
	log.info(invname)
	cmd=" ".join(['$INVFUSION/makeref.real.sh ',invname,bedfile,genbank,fastaref,size,tmp,output])
	log.debug(cmd)
	p=subprocess.call(cmd, shell=True,stderr=subprocess.STDOUT,stdout=file(tmp+"/log."+invname,'w'))
	log.info(p)
	return(p)

def createGenome(log,param):
	f = open(param['list'],'r')
	state=0
	for line in f:
		cols = line.strip().split("\t")
		if (cols[1]=="insilico"):
			e=insilico(cols[0],param['position'],param['href'],param['tempdir'],param['refdir'],param['size'],log)
		else:
			e=real(cols[0],param['position'],cols[2],param['href'],param['tempdir'],param['refdir'],param['size'],log)
		if (e>state):
			state=e
	if (state==0):
		e=mask(param['href'],param['refdir']+"/mask.bed",param['refdir'],param,log)	
	else:
		"error in create virtual genome, stop"
		e=2
	return(e)

def mask(fastaref,maskbed,output,param,log):
	log.info("masking genome")
	cmd=" ".join(['$INVFUSION/mask.genome.sh ',fastaref,maskbed,output])
	p=subprocess.call(cmd, shell=True,stderr=subprocess.STDOUT,stdout=file(param['tempdir']+"/log.mask",'w'))
	return(p)

