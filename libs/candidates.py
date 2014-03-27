
def evalConditions(ts,te,rs,re,con):
	overlap = False
	#print(' '.join([ts,te,rs,re]))
	for c in con:
		#print("##")
		cols = c.split("\t")
		b = eval(cols[1])
		if b:
			return(cols[0])
	return(False)

def applyCondition(inv,con,data,log):
	gtf = data['out'] + "/gtf/" + inv + ".gtf"
	posfile = data['refdir'] + "/" +  inv + ".bed"	
	pos = open(posfile,'r').read().strip().split("\t")
	listkeep = {}
	for line in open(gtf,'r'):
		cols = line.strip().split("\t")
		#log.info(line)
		if cols[2]=="transcript":
			keep = evalConditions(cols[3],cols[4],pos[1],pos[2],con)
			if keep:
				ids = cols[8].split(";")
				listkeep[ids[1]] = keep
	out = open(data['out'] + "/gtf/" + inv + "_filter.gtf",'w')	
	for line in open(gtf,'r'):
		cols = line.strip().split("\t")
		ids = cols[8].split(";")
		if listkeep.has_key(ids[1]):
			tr = ids[1].split(" ")[2].replace('"','')
			exonid = tr + "_" + ids[2].split(" ")[2].replace('"','') 
			out.write(line.strip() + ";" + ' exon_id "' + exonid + '"' +";" + 'condition "' + listkeep[ids[1]] + '"\n' )
	out.close()
	return(0)

def selectCandidates(data,log):
	p = 0 
	f = open(data['conditions'],'r')
	conditions = f.read().strip().split("\n")
	f.close()
	f = open(data['list'],'r')
	for line in f:
		cols = line.strip().split("\t")
		log.info(cols[0])
		p = applyCondition(cols[0],conditions,data,log)
	f.close()
	return(p)

