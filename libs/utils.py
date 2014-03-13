import yaml
import logging


def readconfig(f):
	f = open(f,'r')
	dataMap = yaml.safe_load(f)
	f.close()
	return dataMap
