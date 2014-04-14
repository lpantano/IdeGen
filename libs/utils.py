import yaml
import logging
import subprocess

class showerror( Exception ): pass

def readconfig(f):
	f = open(f,'r')
	dataMap = yaml.safe_load(f)
	f.close()
	return dataMap

def run(name,cmd,err="",out=""):
	if (out!=""):
		p=subprocess.call(cmd, shell=True,stderr=file(err,'a'),stdout=file(out,'a'))
	else:
		p=subprocess.call(cmd, shell=True,stderr=file(err,'a'))
	
	if p!=0:
		raise showerror(name+" broke. See "+err)
	return(p)