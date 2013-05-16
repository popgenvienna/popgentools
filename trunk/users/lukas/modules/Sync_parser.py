import sys
import collections
import numpy as np

#need sync parser
class SyncLineParser(object):
	def __init__(self,syncline, replicates=None, includedel=False):
		'''
		SyncLineParser
		pass it a single line from a sync file
		optional options:
			which (default=None) will indicate which columns have the phenotype of interest
			includedel  (default=False) set to True indicates that the last two columns of sync should be included
		other attributes set on init:
				self.whichcols #
				self.chr 
				self.pos 
				self.ref 
				self.seqs #sync columns with bp information
				self.cmhp
		functions:
			self.get_two_major_alleles()
			 	sets self.majminor, gets the two major alleles, ignoring anything higher than a di-allellic site
			self.get_pop_allele_freqs(two=True)
				gets overall allele freqs from seq information across all replicates
				if two is set to True, major and minor alleles only
			self.get_overall_allele_freqs(two=True)
				gets overall allele freqs from seq information across all replicates
				if two is set to True, major and minor alleles only
			self.get_reduced_allele_freq_dict(two=True)
				gets overall allele freqs from seq information across all replicates
				if two is set to True, major and minor alleles only
			
		'''
		self.includedel = includedel
		if self.includedel:
			self.ncol = 6
		else:
			self.ncol = 4
		self.whichcols=replicates #labels for replicates
		#include last two columns if True
		#parse syncline
		sline = syncline.split()
		self.chr =sline.pop(0)
		self.pos =sline.pop(0)
		self.ref =sline.pop(0)
		self.seqs = [ np.array([int (x) for x in y.split(':')[:self.ncol] ]) for y in sline if ':' in y]
		if ':' not in sline[-1]: #CMH output
			self.cmhp =sline.pop()
		else:
			self.cmhp =None
		#make dictionary with information for phenotype or replicate
		self.majminor = None
		self.pop_allele_freqs =None
		self.overall_allele_freqs=None
		self.seq_dict =None
		self.reduced_seq = None
		self.reduced_seq_dict = None
		self.reduced_af_dict = None


	def get_seq_dict(self):
		if not self.seq_dict:
			self.seq_dict = collections.defaultdict(list)
			for r in range(0, len(self.seqs)):
				self.seq_dict[self.whichcols[r]].append(self.seqs[r])
		return self.seq_dict

	def get_two_major_alleles(self):
		if not self.majminor:
			whichcols = range(0, self.ncol)
			allele_totals = np.array([sum([y[x] for y in self.seqs]) for x in whichcols])
                        # get the highest ranked columns in reverse order
			self.majminor =list(allele_totals.argsort()[-1:-3:-1])
		return self.majminor

	def get_pop_allele_freqs(self, two = True):            
            if not self.pop_allele_freqs:
		if two and not self.majminor:
			self.get_two_major_alleles()
		if two:
			whichcols = self.majminor
		else:
			whichcols = range(0, self.ncol)
                        #print whichcols
		reduced_seq = np.array([[float(x[y]) for y in whichcols] for x in self.seqs]) #reduce
		pop_totals =  [ x.sum() for x in reduced_seq]
		self.pop_allele_freqs =[]
		for i in range(len(self.seqs)):
                    freq="NAN"
                    if (pop_totals[i] > 0):
                        freq=reduced_seq[i]/pop_totals[i]
		    self.pop_allele_freqs.append(freq)
            return self.pop_allele_freqs
	
	def get_overall_allele_freqs(self, two=True):	
		if not self.overall_allele_freqs:
                    self.overall_allele_freqs=[]
		    if not self.pop_allele_freqs:
                        self.get_pop_allele_freqs(two)
                    num_pop=len(self.pop_allele_freqs)
                    self.overall_allele_freqs = [sum([y[x] for y in self.pop_allele_freqs])/num_pop for x in range(0, len(self.pop_allele_freqs[0]))]
		return self.overall_allele_freqs
	
	def get_reduced_seq (self):
		if not self.reduced_seq:
			if not self.majminor:
				self.get_two_major_alleles()
			self.reduced_seq = [[x[y] for y in self.majminor] for x in self.seqs] #reduce 
		return self.reduced_seq

	def get_reduced_seq_dict(self):
		if not self.reduced_seq_dict:
			if not self.reduced_seq:
				self.get_reduced_seq()
			self.reduced_seq_dict = {}
			for r,seq in enumerate(self.reduced_seq):
				rep = self.whichcols[r]
				self.reduced_seq_dict.setdefault(rep,[]).append(seq)
		return self.reduced_seq_dict
            
	def get_reduced_allele_freq_dict(self,two=True):
		if not self.reduced_af_dict:
			if not self.pop_allele_freqs:
				self.get_pop_allele_freqs(two)
			self.reduced_af_dict = collections.defaultdict(list)
			for r,af in enumerate(self.pop_allele_freqs):
				rep = self.whichcols[r]
				self.reduced_af_dict[rep].append(af)
		return self.reduced_af_dict

if __name__ == "__main__":
	import doctest
	doctest.testmod(verbose=1)


