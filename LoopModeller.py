#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import glob
import os
import re

from PDBindexer import PDBindexer

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
		self.ndxedstrands = None
		self.makendxedstrands()
		self.exstrands = None
		self.extendBetaStrands3()
		self.seqstrands = self.deepcopystrands()
		self.extendBetaStrands()
		self.printstrands()
		self.AlignmentFile = self.buildAlignmentFile()
		self.make_template()
		self.loops = None
		self.makeloops()
		
		# executes LoopModelling pipeline with Modeller
		# self.modelBetaBarrel()
		self.modelLoops()
		self.selectModel()

	def deepcopystrands(self):
		tmp = []
		for s in self.strands:
			tmp.append((s[0], s[1]))
		return tmp

	
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
				# the original strand i,j are included and the sequence starts
				# at 1. the -1 and +0 transform them to python array slices
				strands.append((int(strand[0])-1,int(strand[1])+0))
		# print(strands)
		return strands

	def make_template(self):
		pdb = PDBindexer(self.ExtendedStrandsFile)
		resis = pdb.resis
		resicount = 0
		atomcount = 0
		with open(self.TemplateFile, "w") as f:
			for i in range(len(resis)):
				if self.checkifdumping(i):
					resi = resis[i]
					for line in resi:
						self.dump(f, line, atomcount, resicount)
						atomcount += 1
					resicount += 1

	def checkifdumping(self, i):
		dumping = False
		for s in self.exstrands:
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
				l[11:22] +
				"{0: 4d}".format(r+1) +
				l[26:])

	def makendxedstrands(self):
		self.ndxedstrands = []
		n = len(self.strands)
		actualjnexti = 0
		for k in range(n):
			(i, j) = self.strands[k]
			l = j-i + 8
			self.ndxedstrands.append((actualjnexti, actualjnexti+l))
			actualjnexti += l
		# print(self.ndxedstrands)

	def extendBetaStrands3(self):
		n = len(self.strands)
		nextoffset = 0
		seqlen = len(self.sequence)
		self.exstrands = []
		for k in range(n):
			(oi, oj) = self.strands[k]
			(ni, nj) = self.ndxedstrands[k]
			(i, j) = (None, None)
			if k == 0:
				if oi >= 4:
					i = ni
				else:
					i = ni + 4-oi
			else:
				i = ni+4 - nextoffset
			if k == n-1:
				if oj <= seqlen-4:
					j = nj
				else:
					j = nj-seqlen + oj - 4
			else:
				(nextoi, nextoj) = self.strands[k+1]
				lloop = nextoi-oj
				if lloop > 10:
					j = nj
					nextoffset = 4
				elif lloop > 3 and lloop % 2 != 0:
					extension = (lloop-3)/2
					j = nj-4 + extension
					nextoffset = extension
				elif lloop > 2 and lloop % 2 == 0:
					extension = (lloop-2)/2
					j = nj-4 + extension
					nextoffset = extension
				else:
					j = nj-4
					nextoffset = 0
			self.exstrands.append((i, j))
		# print(self.exstrands)

	def printstrands(self):
		for i in range(len(self.strands)):
			print(self.strands[i])
			print(self.seqstrands[i])
			print(self.ndxedstrands[i])
			print(self.exstrands[i])
			print("")
		print(len(self.sequence))

	
	def extendBetaStrands(self):
		n = len(self.seqstrands)
		for i in range(n): # access by index to prevent weird behaviour when modifying next strands
			(beg_cur, end_cur) = self.seqstrands[i]
			if i == 0:
				if beg_cur >= 4:
					beg_cur -= 4
				else:
					beg_cur = 0
			if i < n-1:
				(beg_nex, end_nex) = self.seqstrands[i+1]
				dist = beg_nex - end_cur

				# ß-barrel extension
				# loop longer than 10 resideus
				if dist > 10:
					end_cur = end_cur + 4
					beg_nex = beg_nex - 4
					dist = dist - 8
					self.seqstrands[i] = (beg_cur, end_cur)
					self.seqstrands[i+1] = (beg_nex, end_nex)
				# odd-loop shorter than 10 but longer than 2
				elif dist > 3 and dist % 2 != 0:
					# odd numbered loop
					while dist > 3:
						end_cur = end_cur + 1
						beg_nex = beg_nex - 1
						dist = dist - 2

					self.seqstrands[i] = (beg_cur, end_cur)
					self.seqstrands[i+1] = (beg_nex, end_nex)
				
				elif dist > 2 and dist % 2 == 0:
					# even numbered loop	
					while dist > 2:
						end_cur = end_cur + 1
						beg_nex = beg_nex - 1
						dist = dist - 2

					self.seqstrands[i] = (beg_cur, end_cur)
					self.seqstrands[i+1] = (beg_nex, end_nex)			
				
				else:
					pass
			else:
				if end_cur < n:
					if end_cur <= n-4:
						end_cur += 4
					else:
						end_cur = n
					self.seqstrands[i] = (beg_cur, end_cur)
		#print self.seqstrands



	
	def buildAlignmentFile(self):
		# fills alignment[] with gaps
		alignment = []
		for seq in self.sequence:
			alignment.append("-")

		# fill the structure with sequence when stranded (according to self.strands )
		for (strand_beg, strand_end) in self.seqstrands:
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

	def makeloops(self):
		n = len(self.seqstrands)
		loops = []
		for k in range(n):
			(i, j) = self.seqstrands[k]
			if k == 0:
				if i > 0:
					loops.append((1, i))
			elif k == n-1:
				if j < len(self.sequence)-1:
					loops.append((j+1, len(self.sequence)))
			else:
				nexti = self.seqstrands[k+1][0]
				loops.append((i, j+1))
		self.loops = loops

	def formatLoopResidues(self):
		formattedRes = ""
		n = len(self.loops)
		for k in range(n):
			(i, j) = self.loops[k]
			s = "self.residue_range('{0}:A','{1}:A')".format(i,j)
			formattedRes += s
			if k < n-1:
				formattedRes += ','
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
		# aBetaBerrelModel.md_level = refine.very_slow
		aBetaBerrelModel.md_level = refine.very_fast
		# model ß-barrel
		aBetaBerrelModel.make()


	def modelLoops(self):
		selectedResidues = self.formatLoopResidues()	
		print(selectedResidues)
		# Create a new class based on 'loopmodel' so that we can redefine
		# select_loop_atoms (necessary)
		class MyLoop(loopmodel):
		    # This routine picks the residues to be refined by loop modeling
		    def select_loop_atoms(self):
		        return selection( eval(selectedResidues) )
			def special_patches(self, aln):
					for chain in self.chains:
						chain.name = 'A'		
		# define ß-berrel loops' modelling parameters
		# aLoop = MyLoop(self.ModellerEnv, inimodel=self.TemplateFile, sequence=self.FastaFile, loop_assess_methods=assess.DOPE)
		aLoop = MyLoop(self.ModellerEnv, inimodel="A9WGN5_full.B99990001.pdb", sequence=self.FastaFile, loop_assess_methods=assess.DOPE)
		aLoop.loop.starting_model = 1
		aLoop.loop.ending_model = self.nModels
		# define loops modelling refinement level
		# aLoop.loop.md_level = refine.very_slow
		aLoop.loop.md_level = refine.very_fast
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

