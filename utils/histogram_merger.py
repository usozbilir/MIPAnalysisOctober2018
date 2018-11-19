import os
from copy import deepcopy
import ROOT
import time

#class merging ROOT files containing addable objects within different files
#usage:
#        from histogram_merger import HistogramMerger
#        with HistogramMerger(output_path) as hm:
#            hm.addHistogramFile(input1.path)
#            hm.addHistogramFile(input2.path)
#            hm.addHistogramFile(input3.path)
                



class HistogramMerger(object):
	
	def __init__(self, path):
		super(HistogramMerger, self).__init__()
		self.path = path
		self.sources = []
		self.merged = False
		self.objects = []
		self.objectsToBeMerged = []

		self.mergedObjects = {}
		
		self.patterns = []

		self.currentROOTFile = None


	def __enter__(self):
		return self

	def __exit__(self ,type, value, traceback):
		if self.merged==True:
			return

		self.merged=True
		self.mergeHistograms()
		self.closeCurrentROOTFile()
		return False

	def addSelection(self, regexp):	#object must fulfill the indicated pattern as given as regexp.
		self.patterns.append(regexp)

	def addHistogramFile(self, sourcePath):
		if not os.path.exists(sourcePath):
			print "%s does not exist and is not added" % sourcePath
			return			
		if sourcePath in self.sources:
			print "%s has already been added" % sourcePath
			return

		self.sources.append(sourcePath)
		
		if len(self.sources)==1:
			directories = {}
			self.loadObjects(directories, sourcePath)
			print "%s objects in first file." % len(self.objects)
			
	def mergeHistograms(self):
		if len(self.patterns) > 0:
			import re
			for name in self.objects:
				for pattern in self.patterns:
					if bool(re.match(pattern, name)):
						self.objectsToBeMerged.append(name)
						break
		else:
			self.objectsToBeMerged = self.objects

		self.objectsToBeMerged = sorted(self.objectsToBeMerged)

		Ns = len(self.sources)
		No = len(self.objectsToBeMerged)
		print "%s objects to be merged across %s files." % (No, Ns) 
		
		startTime = time.time()
		self.mergedObjects = {}
		for s in range(Ns):
			source = self.sources[s]
			self.loadHistogramsFromCurrentFile(source)

			print "After %s seconds: %s/%s files merged" % (time.time()-startTime, s+1, Ns)

		self.openROOTFile(self.path, "RECREATE")
		for name in self.mergedObjects:
			filedir = "/".join(name.split("/")[0:-1])
			if not self.currentROOTFile.GetDirectory(filedir):
				self.currentROOTFile.mkdir(filedir)
			self.currentROOTFile.cd(filedir)
			self.mergedObjects[name].Write()
		self.closeCurrentROOTFile()


	def loadHistogramsFromCurrentFile(self, source):
		self.openROOTFile(source)
		N = len(self.objectsToBeMerged)
		for i in range(N):
			name = self.objectsToBeMerged[i]
			loaded_histogram = deepcopy(self.currentROOTFile.Get(name))
			if not hasattr(loaded_histogram, "Add"):
				print name,"is not addable."
				continue
			if name not in self.mergedObjects:
				self.mergedObjects[name] = loaded_histogram
			else:
				try:
					self.mergedObjects[name].Add(loaded_histogram)
				except:
					print name,"is not mergeable."
					continue


	def loadObjects(self, dir_keys, filePath=None):
		if not filePath==None:
			self.openROOTFile(filePath)			
			for key in self.currentROOTFile.GetListOfKeys():
				dir_keys[key.GetName()] = {}
			self.loadObjects(dir_keys)
		else:
			for key in dir_keys:
				self.currentROOTFile.cd()
				if self.currentROOTFile.GetDirectory(key):
					self.currentROOTFile.cd(key)
				else:
					self.objects.append(key)
					return
				for new_key in ROOT.gDirectory.GetListOfKeys(): 
					if self.currentROOTFile.GetDirectory("%s/%s"%(key, new_key.GetName())):
						dir_keys[key]["%s/%s"%(key, new_key.GetName())] = {}
					else: 
						self.objects.append("%s/%s"%(key, new_key.GetName()))
							
				self.loadObjects(dir_keys[key])


	def openROOTFile(self, filePath, option="READ"):
		self.closeCurrentROOTFile()
		self.currentROOTFile = ROOT.TFile(filePath, option)

	def closeCurrentROOTFile(self):
		if not self.currentROOTFile==None:
			self.currentROOTFile.Close()
			self.currentROOTFile = None
