#!/usr/bin/env python
# -*- coding: utf-8 -*-
# import subprocess
import argparse
# import shlex
import glob
import os.path
import re

from modeller import *
from modeller.automodel import *


class LoopModeller:

	def __init__(self, FastaID, FastaFile, StrandsFile, TemplateFile, isScrambled, nModels):
		# define class attributes
		self.FastaID = FastaID
		self.FastaFile = FastaFile
		self.StrandsFile = StrandsFile
		self.TemplateFile = TemplateFile
		self.isScrambled = isScrambled
		self.nModels = nModels
		self.BasePath = os.getcwd()
		self.ModellerEnv = environ()
		self.sequence = self.readSequenceFromFASTA()
		self.numbering = self.readPDBnumbering()
		self.strands = self.readStrandsFile()
		self.extendBetaStrands()
		self.AlignmentFile = self.buildAlignmentFile()
		
		# executes LoopModelling pipeline with Modeller
		self.modelBetaBarrel()
		self.modelLoops()
		self.selectModel()

	
	def readSequenceFromFASTA(self):
		seq = []
		with open(self.FastaFile) as f:
			lines = f.readlines()
			for l in lines:
				if not re.search(r'^>', l):
					chars = re.search(r'([A-Za-z]+)', l).group(0)
					seq.extend(chars)
		# print(seq)
		return seq

	
	def readStrandsFile(self):
		strands = []
		with open(self.StrandsFile) as f:
			lines = f.readlines()
			for l in lines:
				strand = re.search(r'(\d+)\s+(\d+)',l).groups(0)
				strands.append(strand)
		# print(strands)
		return strands

	
	def readPDBnumbering(self):
		numbering = []
		with open(self.TemplateFile) as f:
			lines = f.readlines()
			for l in lines:
				if re.search(r'^ATOM', l) and re.search(r'CA', l):
					number = int(l[22:26])
					numbering.append(number)
					for it in iter(self.sequence[1:]):
						number = number + 1 
						numbering.append(number)
					break
		# print(numbering)
		return numbering
					

	
	def extendBetaStrands(self):
		for i, (beg_cur, end_cur) in enumerate(self.strands):
			if i < len(self.strands)-1:
				(beg_nex, end_nex) = self.strands[i+1]
				dist = int(beg_nex)-int(end_cur) - 1 

				# ß-barrel extension
				# loop longer than 10 resideus
				if dist > 10:
					end_cur = int(end_cur) + 4
					beg_nex = int(beg_nex) - 4
					self.strands[i] = (beg_cur, end_cur)
					self.strands[i+1] = (beg_nex, end_nex)
				# odd-loop shorter than 10 but longer than 2
				elif dist > 3 and dist % 2 > 0:
					# odd numbered loop
					while dist > 3:
						end_cur = int(end_cur) + 1
						beg_nex = int(beg_nex) - 1
						dist = dist - 2

					self.strands[i] = (beg_cur, end_cur)
					self.strands[i+1] = (beg_nex, end_nex)
				
				elif dist > 2 and dist % 2 == 0:
					# even numbered loop	
					while dist > 2:
						end_cur = int(end_cur) + 1
						beg_nex = int(beg_nex) - 1
						dist = dist - 2

					self.strands[i] = (beg_cur, end_cur)
					self.strands[i+1] = (beg_nex, end_nex)			
				
				else:
					pass
		# print self.strands

	
	def buildAlignmentFile(self):
		# fills alignment[] with gaps
		alignment = []
		for seq in self.sequence:
			alignment.append("-")
		# fill the structure with sequence when stranded (according to self.strands )
		for (strand_beg, strand_end) in self.strands:
			for i in range( (int(strand_beg) - 1 - self.numbering[0]), (int(strand_end) - self.numbering[0]) ):
				alignment[i] = self.sequence[i]
		# output alignment *.pir file
		AlignmentFile = str(self.FastaID) + ".pir"
		with open(AlignmentFile, "w") as f:
			f.write( ">P1;{0}\n".format(self.FastaID) )
			f.write( "structure:{0}:{1}:A: : : : : :\n".format(self.TemplateFile, self.numbering[0]) )
			f.write( "{0}*\n".format(''.join(self.sequence)) )
			f.write( "\n" )

			f.write( ">P1;{0}_full\n".format(self.FastaID) )
			f.write( "structure:{0}:{1}:A: : : : : :\n".format(self.TemplateFile, self.numbering[0]) )
			f.write( "{0}*\n".format(''.join(alignment)) )
		
		return AlignmentFile


	def assignChain(self, chain):
		pass

	def modelBetaBarrel(self):
		pass

	def modelLoops(self):
		pass

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
		elif not re.search(r'(\S+)\.fa', FastaID):
			raise ValueError("The FASTA files must contain a FASTA identifier ID")
		else:
			FastaID = re.search(r'^(\S+).*\.fa', FastaID).group(1)
		
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
		# LoopModelScrambled = LoopModeller(FastaID, ScrambledFastaFile, StrandsFile, TemplateFile, True, nModels)
		# LoopModels.append( (LoopModel,LoopModelScrambled) )

