#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import glob
import os
import re
# from Bio.PDB import *
from modeller import *
from modeller.automodel import *


class LoopModeller:

	def __init__(self, FastaID, FastaFile, StrandsFile, TemplateFile, isScrambled, nModels):
		# define class attributes
		self.FastaID = FastaID
		self.FastaFile = FastaFile
		self.StrandsFile = StrandsFile
		self.ExtendedStrandsFile = TemplateFile
		self.TemplateFile = "{0}_renumbered_extended.pdb".format(self.FastaID)
		self.isScrambled = isScrambled
		self.nModels = nModels
		self.BasePath = os.getcwd()
		self.ModellerEnv = environ()
		self.numbering = [] # self.numbering will be filled form [1,N] in readSequenceFromFASTA
		self.sequence = self.readSequenceFromFASTA()
		self.strands = self.readStrandsFile()
		self.extendBetaStrands()
		self.AlignmentFile = self.buildAlignmentFile()
		self.make_template()
		# self.resi = self.readPDBnumbering()
		
		# executes LoopModelling pipeline with Modeller
		self.modelBetaBarrel()
		self.modelLoops()
		self.selectModel()

	
	def readSequenceFromFASTA(self):
		seq = []
		firstTagSeen = False
		secondTagSeen = False
		with open(self.FastaFile) as f:
			lines = f.readlines()
			for l in lines:
				if re.search(r'^>',l) and not firstTagSeen:
					firstTagSeen = True
				elif not re.search(r'^>(tr|)*(\S+)', l) and not secondTagSeen:
					chars = re.search(r'([A-Za-z]+)', l).group(0)
					seq.extend(chars)
				elif re.search(r'^>',l) and firstTagSeen and not secondTagSeen:
					secondTagSeen = True
					break
		# build corresbonding self.numbering sequence
		i = 1
		for c in seq:
			self.numbering.append(i)
			i = i + 1
		# print(seq)
		return seq

	
	def readStrandsFile(self):
		strands = []
		with open(self.StrandsFile) as f:
			lines = f.readlines()
			for l in lines:
				strand = re.search(r'(\d+)\s+(\d+)',l).groups(0)
				strands.append((int(strand[0]),int(strand[1])))
		print(strands)
		return strands

	
	# def readPDBnumbering(self):
	# 	resi = []
	# 	with open(self.TemplateFile) as f:
	# 		lines = f.readlines()
	# 		for i,l in enumerate(lines):
	# 			if re.search(r'^ATOM', l) and re.search(r'CA', l):
	# 				resnum = int(l[22:26])
	# 				resnam = "{0:3}".format(l[17:20])
	# 				chain  = l[21]
	# 				resi.append( (resnam, resnum, chain, l) )
	# 	for (bs, es) in self.strands:
	# 		strandSeq = "".join(self.sequence[int(bs):int(es)])
	# 		for i, (res, num, cha, lin) in enumerate(resi):
				
				
	# 	return resi
					
	def make_template(self):
		with open(self.ExtendedStrandsFile) as f:
			lines = f.readlines()
		with open(self.TemplateFile, "w") as f:
			atomcount = 0
			actresicount = 0
			resicount = 0
			lastresnum = None
			first = True
			dumping = False
			n = len(lines)
			for i in range(n):
				l = lines[i]
				nl = None
				if i < n-1:
					nl = lines[i+1]
				if (l[0:4] == 'ATOM' or l[0:4] == 'HETA'):
					resnum = int(l[22:26])
					print(resnum)
					nextresnum = None
					if (nl[0:4] == 'ATOM' or nl[0:4] == 'HETA'):
						nextresnum = int(nl[22:26])
						print(nextresnum)
					if first:
						print("here1")
						dumping = self.checkifdumping(actresicount)
						if dumping:
							self.dump(f, l, atomcount, resicount)
							atomcount += 1
							first = False
					elif (resnum == nextresnum):
						print("here2")
						if dumping:
							self.dump(f, l, atomcount, resicount)
							atomcount += 1
					elif nextresnum is not None:
						print("here3")
						dumping = self.checkifdumping(actresicount)
						if dumping:
							self.dump(f, l, atomcount, resicount)
							atomcount += 1
							resicount += 1
						actresicount += 1
					else:
						print("here4")
						dumping = self.checkifdumping(actresicount)
						if dumping:
							self.dump(f, l, atomcount, resicount)
					lastresnum = resnum

	def checkifdumping(self, i):
		dumping = False
		for s in self.strands:
			if i >= s[0] and i < s[1]:
				dumping = True
		return(dumping)
					
	def dump(self, f, l, a, r):
		"""
		f is filehandle
		l original line to dump
		a is current atom number (starting at 0)
		r is current resi number (starting at 0)
		"""
		f.write(l[0:6] +
				"{0: 5d}".format(a+1) +
				l[12:22] +
				"{0: 4d}".format(r+1) +
				l[26:])

	
	def extendBetaStrands(self):
		n = len(self.strands)
		for i in range(n): # access by index to prevent weird behaviour when modifying next strands
			(beg_cur, end_cur) = self.strands[i]
			if i == 0:
				if beg_cur >= 4:
					beg_cur -= 4
				else:
					beg_cur = 0
			if i < n-1:
				(beg_nex, end_nex) = self.strands[i+1]
				dist = beg_nex - end_cur

				# ß-barrel extension
				# loop longer than 10 resideus
				if dist > 10:
					end_cur = end_cur + 4
					beg_nex = beg_nex - 4
					dist = dist - 8
					self.strands[i] = (beg_cur, end_cur)
					self.strands[i+1] = (beg_nex, end_nex)
				# odd-loop shorter than 10 but longer than 2
				elif dist > 3 and dist % 2 != 0:
					# odd numbered loop
					while dist > 3:
						end_cur = end_cur + 1
						beg_nex = beg_nex - 1
						dist = dist - 2

					self.strands[i] = (beg_cur, end_cur)
					self.strands[i+1] = (beg_nex, end_nex)
				
				elif dist > 2 and dist % 2 == 0:
					# even numbered loop	
					while dist > 2:
						end_cur = end_cur + 1
						beg_nex = beg_nex - 1
						dist = dist - 2

					self.strands[i] = (beg_cur, end_cur)
					self.strands[i+1] = (beg_nex, end_nex)			
				
				else:
					pass
			else:
				if end_cur < n:
					if end_cur <= n-4:
						end_cur += 4
					else:
						end_cur = n
					self.strands[i] = (beg_cur, end_cur)
		print self.strands

	
	def buildAlignmentFile(self):
		# fills alignment[] with gaps
		alignment = []
		for seq in self.sequence:
			alignment.append("-")

		# fill the structure with sequence when stranded (according to self.strands )
		for (strand_beg, strand_end) in self.strands:
			for i in range(strand_beg, strand_end):
				alignment[i] = self.sequence[i]
		# output alignment *.pir file
		AlignmentFile = "{0}/{1}.pir".format(self.BasePath, self.FastaID)
		with open(AlignmentFile, "w") as f:
			f.write( ">P1;{0}\n".format(self.FastaID) )
			f.write( "structure:{0}:{1}:A: : : : : :\n".format(self.TemplateFile, self.numbering[0]) )
			f.write( "{0}*\n".format(''.join(alignment)) )
			f.write( "\n" )

			f.write( ">P1;{0}_full\n".format(self.FastaID) )
			f.write( "sequence:{0}:{1}:A: : : : : :\n".format(self.FastaID+"_full", self.numbering[0]) )
			f.write( "{0}*\n".format(''.join(self.sequence)) )
		if os.path.exists(AlignmentFile):
			return AlignmentFile
		else:
			raise ValueError("Unable to create modellers' alignment (*.pir) file")



	def formatLoopResidues(self):
		formattedRes = ""
		memory_value = 0
		# build formatted strings of from self.strands
		for(beg_strand, end_strand) in self.strands:
			if memory_value == 0:
				s = "self.residue_range('1:A','{0}:A')".format(beg_strand)
				memory_value = int(end_strand)
			else:
				s = "self.residue_range('{0}:A','{1}:A')".format(memory_value, int(beg_strand))
				memory_value = int(end_strand)
			formattedRes.join(s)
		# join the last loop of the model to the selection
		formattedRes.join("self.residue_range('{0}:A','{1}:A')".format(memory_value, len(self.sequence)))
		print formattedRes
		# return the string built from selected (all) loops
		return formattedRes


	def modelBetaBarrel(self):
		# redefine modeller's model class to rename chain
		class MyModel(automodel):
			def special_patches(self, aln):
				for chain in self.chains:
					chain.name = 'A'

		# set ModellerEnvironment parameters
		self.ModellerEnv.io.atom_files_directory = [self.BasePath]
		self.ModellerEnv.schedule_scale = physical.values(default=1.0, soft_sphere=0.7)
		# define ß-berrel modelling parameters
		aBetaBerrelModel = MyModel(self.ModellerEnv, alnfile=self.AlignmentFile, knowns=self.FastaID, sequence=self.FastaID + "_full")
		aBetaBerrelModel.starting_model = 1                 # index of the first model
		aBetaBerrelModel.ending_model = self.nModels            # index of the last model
		# simulation parameters
		aBetaBerrelModel.library_schedule = autosched.slow
		aBetaBerrelModel.max_var_iterations = 300
		aBetaBerrelModel.repeat_optimization = 20
		aBetaBerrelModel.max_molpdf = 1e6
		# define refinement level
		aBetaBerrelModel.md_level = refine.very_slow
		# model ß-barrel
		aBetaBerrelModel.make()


	def modelLoops(self):
		selectedResidues = self.formatLoopResidues()	
		
		# Create a new class based on 'loopmodel' so that we can redefine
		# select_loop_atoms (necessary)
		class MyLoop(loopmodel):
		    # This routine picks the residues to be refined by loop modeling
		    def select_loop_atoms(self):
		        return selection(selectedResidues)
			def special_patches(self, aln):
					for chain in self.chains:
						chain.name = 'A'		
		# define ß-berrel loops' modelling parameters
		aLoop = MyLoop(self.ModellerEnvironment, inimodel=self.TemplateFile, sequence=self.FastaFile, loop_assess_methods=assess.DOPE)
		aLoop.loop.starting_model = 1
		aLoop.loop.ending_model = self.nModels
		# define loops modelling refinement level
		aLoop.loop.md_level = refine.very_slow
		# model loops
		aLoop.make()

	def selectModel(self):
		pass

#####################################################################
# define command line arguments
ArgParser = argparse.ArgumentParser()
ArgParser.add_argument("-f", "--fasta", "--fasta-file", "--fastafile", metavar = "FASTA", type = str, help="Path to the <FASTA file>", required=True)
ArgParser.add_argument("-s", "--strands", "--strands-file", "--strandsfile", metavar = "STRANDS", type = str, help="Path to the <STRANDS file>", required=True)
ArgParser.add_argument("-t", "--template", "--PDB", "--template-file", "--PDB-file", metavar = "PDB", type=str, required=True, help="Path to the <PDB template file>")
ArgParser.add_argument("-n", "--nModels", "--nModels", metavar="N", type=int, required=True, help="Number of models to be generated")
ArgParser.add_argument("-v", "--verbose", action="store_true", help="Activates Modeller's verbose mode")

# main execution of this script
if __name__ == "__main__":
	# parse command line arguments
	args = ArgParser.parse_args()
	
	FastaFiles = []
	file = args.fasta
	if re.search(r'\*', file):
		FastaFiles = glob.glob(file)
		StrandsFiles = glob.glob("*.strands")
		PDBFiles = glob.glob("*_ext_l04.pdb")
	if len(FastaFiles) < 1:
		if os.path.exists(file):
			FastaFiles.append(file)
		else:
			raise ValueError("File {} does not exists.".format(FastaFile))

	LoopModels = []
	# loop each FastaFile
	for FastaFile in FastaFiles:
	
		# specifies *.fa or *.fasta file
		FastaID = os.path.split(FastaFile)[1]
		if not isinstance(FastaID, str):
			raise TypeError("The FASTA ID ({}) must be a string".format(FastaID))
		elif not re.search(r'(\S+)\.pfa', FastaID) and not re.search(r'(\S+)\.fa', FastaID):
			raise ValueError("The FASTA files must contain a FASTA identifier ID")
		else:
			FastaID = FastaID[0:6]

		# look form scrambled FASTA file and TEMPLATE PDB file
		# print(FastaFile)
		# ScrambledFastaFile = re.sub(FastaID,"{0}_scrambled".format(FastaID), os.path.split(FastaFile)[1])
		# print(ScrambleFastaFile)
		# if not os.path.exists(ScrambledFastaFile):
		# 	raise ValueError("Unable to read the scrambled FASTA file {0}".format(ScrambledFastaFile))
		
		# ScrambledTemplateFile = 
		
		# specifies the *.strands file
		StrandsFile = args.strands
		# StrandsFile = re.sub(r'(\.fa)', r'.strands', FastaFile)
		if not os.path.exists(StrandsFile):
			raise ValueError("File {} does not exists.".format(StrandsFile))
		
		# specifies the PDB template file
		TemplateFile = args.template
		if not os.path.exists(TemplateFile):
			raise ValueError("File {} does not exists.".format(TemplateFile))
		
		# specifies mumber of models to generate
		nModels = args.nModels
		if nModels < 1:
			raise ValueError("The number of requested models ({}) must be at lest one".format(nModels))

		# activates Modeller's verbose mode
		if args.verbose:
			log.verbose()

		LoopModel = LoopModeller(FastaID, FastaFile, StrandsFile, TemplateFile, False, nModels)
		#LoopModelScrambled = LoopModeller(FastaID, ScrambledFastaFile, StrandsFile, ScrambledTemplateFile, True, nModels)
		#LoopModels.append( (LoopModel,LoopModelScrambled) )

