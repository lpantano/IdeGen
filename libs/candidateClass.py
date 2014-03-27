class candidate:
	def __init__(self,idx):
		self.id=idx
		self.tr = "none"
		self.i1reads = 0
		self.i2reads = 0
		self.i1bp = "NO"
		self.i2bp = "NO"
		self.known = 0
		self.cpx = 0
		self.cf = 0
		self.g1 = []
		self.g2 = []
		self.condition = "none"
