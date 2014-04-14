class candidate:
	def __init__(self,idx):
		self.id=idx
		self.number=self.getnumber(idx)
		self.tr = "none"
		self.i1reads = 0
		self.i2reads = 0
		self.i1bp = "NO"
		self.i2bp = "NO"
		self.known = 0
		self.cpx = 0
		self.cf = 0.0
		self.g1 = [0,0,0]
		self.g2 = [0,0,0]
		self.condition = "none"
	def getnumber(self,i):
		return(int(i.split("_")[1]))
	def setcf(self):
		self.cf = ((1.0*self.g1[2]+1)/(1.0*self.g2[1]+1))
		return(self.cf)