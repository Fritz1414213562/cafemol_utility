# compiler
CC=g++
#CFLAGS= -std=c++17 -pedantic -Wall -g -v
CFLAGS= -std=c++17 -pedantic -Wall -g 
INCLUDE= -I. -I../include -I../include/misc

# directory
DCDDIR= ../include/DCD
PDBDIR= ../include/PDB
NINFODIR= ../include/NINFO
PSFDIR= ../include/PSF
MISCDIR= ../include/misc
EXEDIR= ../bin
OBJDIR= ./obj

# generate exe file
#SLIDE= slideAnalyze.o FileOpenJudge.o DCDAnalyzer.o DCDParser.o PDBReader.o PDBParser.o
#DCD2PDB= makePDBfromDCD.o FileOpenJudge.o DCDAnalyzer.o DCDParser.o PDBWriter.o PDBParser.o
#SEARCHNATIVE= searchNativeContact.o FileOpenJudge.o DCDAnalyzer.o DCDParser.o PDBReader.o PDBParser.o
#NCRANGE= nativeContactedRange.o FileOpenJudge.o DCDAnalyzer.o DCDParser.o PDBReader.o PDBParser.o
#HISTOGRAM= convertHistData.o DensityFunction.o
#QSCORE= FileOpenJudge.o NinfoReader.o NinfoParser.o DCDAnalyzer.o DCDParser.o
#SEARCHCONTACTHISTONE= searchContactedHistoneChain.o FileOpenJudge.o DCDAnalyzer.o DCDParser.o PSFReader.o PSFParser.o
#GETPDNS= getProteinResiContacted2DNA.o FileOpenJudge.o NinfoReader.o NinfoParser.o PSFReader.o PSFParser.o PDBReader.o PDBParser.o
#UNWRAPEXTENT= searchUnwrappedExtent.o FileOpenJudge.o DCDAnalyzer.o DCDParser.o PSFReader.o PSFParser.o NinfoReader.o NinfoParser.o
#SLIDE2DYAD= searchDNAContacted2Dyad.o DCDAnalyzer.o DCDParser.o PSFReader.o PSFParser.o
#HELICASEPROCESS= helicaseProcessingExtent.o DCDAnalyzer.o DCDParser.o PSFReader.o PSFParser.o
#PDNSQSCORE= calcPDNSQscore.o NinfoReader.o NinfoParser.o DCDAnalyzer.o DCDParser.o PSFReader.o PSFParser.o
#SORTPDNS= sortPDNSResidue.o DCDAnalyzer.o DCDParser.o PSFReader.o PSFParser.o NinfoReader.o NinfoParser.o
#DETECTDNAPROCONTACT= detectDNAProteinContacts.o DCDAnalyzer.o DCDParser.o PSFReader.o PSFParser.o
#PDNSMAPPING= generateDNAProteinHeatMap.o DCDAnalyzer.o DCDParser.o PSFReader.o PSFParser.o
#CONTACTSTATE2PROB= calcContactStateProbability.o PSFReader.o PSFParser.o
#GETCONTACTADDRESS= getPDNS_ELEContacts.o PSFReader.o PSFParser.o PDBReader.o PDBParser.o
#SHOWSTATES= showStateTransition.o FileOpenJudge.o
#TEST= DCDWriter.o DCDParser.o

SLIDE= $(addprefix $(OBJDIR)/, slideAnalyze.o FileOpenJudge.o DCDAnalyzer.o DCDParser.o PDBReader.o PDBParser.o)
DCD2PDB= $(addprefix $(OBJDIR)/, makePDBfromDCD.o FileOpenJudge.o DCDReader.o DCDParser.o PDBWriter.o PDBParser.o)
SEARCHNATIVE= $(addprefix $(OBJDIR)/, searchNativeContact.o FileOpenJudge.o DCDAnalyzer.o DCDParser.o PDBReader.o PDBParser.o)
NCRANGE= $(addprefix $(OBJDIR)/, nativeContactedRange.o FileOpenJudge.o DCDAnalyzer.o DCDParser.o PDBReader.o PDBParser.o)
HISTOGRAM= $(addprefix $(OBJDIR)/, convertHistData.o DensityFunction.o)
QSCORE= $(addprefix $(OBJDIR)/, calcQscore.o NinfoReader.o NinfoParser.o DCDAnalyzer.o DCDParser.o)
SEARCHCONTACTHISTONE= $(addprefix $(OBJDIR)/, searchContactedHistoneChain.o FileOpenJudge.o DCDAnalyzer.o DCDParser.o PSFReader.o PSFParser.o)
GETPDNS= $(addprefix $(OBJDIR)/, getProteinResiContacted2DNA.o FileOpenJudge.o NinfoReader.o NinfoParser.o PSFReader.o PSFParser.o PDBReader.o PDBParser.o)
UNWRAPEXTENT= $(addprefix $(OBJDIR)/, searchUnwrappedExtent.o FileOpenJudge.o DCDAnalyzer.o DCDParser.o PSFReader.o PSFParser.o NinfoReader.o NinfoParser.o)
SLIDE2DYAD= $(addprefix $(OBJDIR)/, searchDNAContacted2Dyad.o DCDAnalyzer.o DCDParser.o PSFReader.o PSFParser.o)
HELICASEPROCESS= $(addprefix $(OBJDIR)/, helicaseProcessingExtent.o DCDAnalyzer.o DCDParser.o PSFReader.o PSFParser.o)
PDNSQSCORE= $(addprefix $(OBJDIR)/, calcPDNSQscore.o NinfoReader.o NinfoParser.o DCDAnalyzer.o DCDParser.o PSFReader.o PSFParser.o)
SORTPDNS= $(addprefix $(OBJDIR)/, sortPDNSResidue.o DCDAnalyzer.o DCDParser.o PSFReader.o PSFParser.o NinfoReader.o NinfoParser.o)
DETECTDNAPROCONTACT= $(addprefix $(OBJDIR)/, detectDNAProteinContacts.o DCDAnalyzer.o DCDParser.o PSFReader.o PSFParser.o)
PDNSMAPPING= $(addprefix $(OBJDIR)/, generateDNAProteinHeatMap.o DCDAnalyzer.o DCDParser.o PSFReader.o PSFParser.o)
CONTACTSTATE2PROB= $(addprefix $(OBJDIR)/, calcContactStateProbability.o PSFReader.o PSFParser.o)
GETCONTACTADDRESS= $(addprefix $(OBJDIR)/, getPDNS_ELEContacts.o PSFReader.o PSFParser.o PDBReader.o PDBParser.o)
SHOWSTATES= $(addprefix $(OBJDIR)/, showStateTransition.o FileOpenJudge.o)
TEST= $(addprefix $(OBJDIR)/, test.o DCDWriter.o DCDParser.o)


slide: $(SLIDE)
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix $(EXEDIR)/, slide) $(SLIDE)

dcd2pdb: $(DCD2PDB)
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix $(EXEDIR)/, dcd2pdb) $(DCD2PDB)

searchnative: $(SEARCHNATIVE)
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix $(EXEDIR)/, searchnative) $(SEARCHNATIVE)

ncrange: $(NCRANGE)
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix $(EXEDIR)/, ncrange) $(NCRANGE)

histogram: $(HISTOGRAM)
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix $(EXEDIR)/, histogram) $(HISTOGRAM)

qscore: $(QSCORE)
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix $(EXEDIR)/, qscore) $(QSCORE)

searchcontacthistone: $(SEARCHCONTACTHISTONE)
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix $(EXEDIR)/, searchcontacthistone) $(SEARCHCONTACTHISTONE)

getpdns: $(GETPDNS)
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix $(EXEDIR)/, getpdns) $(GETPDNS)

unwrapextent: $(UNWRAPEXTENT)
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix $(EXEDIR)/, unwrapextent) $(UNWRAPEXTENT)

slide2dyad: $(SLIDE2DYAD)
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix $(EXEDIR)/, slide2dyad) $(SLIDE2DYAD)

helicaseprocess: $(HELICASEPROCESS)
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix $(EXEDIR)/, helicaseprocess) $(HELICASEPROCESS)

pdnsqscore: $(PDNSQSCORE)
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix $(EXEDIR)/, pdnsqscore) $(PDNSQSCORE)

sortpdns: $(SORTPDNS)
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix $(EXEDIR)/, sortpdns) $(SORTPDNS)

detectdnaprocontact: $(DETECTDNAPROCONTACT)
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix $(EXEDIR)/, detectdnaprocontact) $(DETECTDNAPROCONTACT)

pdnsmapping: $(PDNSMAPPING)
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix $(EXEDIR)/, pdnsmapping) $(PDNSMAPPING)

contactstate2prob: $(CONTACTSTATE2PROB)
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix $(EXEDIR)/, contactstate2prob) $(CONTACTSTATE2PROB)

getcontactaddress: $(GETCONTACTADDRESS)
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix $(EXEDIR)/, getcontactaddress) $(GETCONTACTADDRESS)

showstates: $(SHOWSTATES)
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix $(EXEDIR)/, showstates) $(SHOWSTATES)

# clean
clean:
	rm obj/*.o

test: $(TEST)
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix $(EXEDIR)/, test) $(TEST)

# object file

$(addprefix $(OBJDIR)/, DCDParser.o): $(addprefix $(DCDDIR)/, DCDParser.cpp)
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix $(OBJDIR)/, DCDParser.o) -c $(addprefix $(DCDDIR)/, DCDParser.cpp) 

$(addprefix $(OBJDIR)/, DCDAnalyzer.o): $(addprefix $(DCDDIR)/, DCDAnalyzer.cpp)
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix $(OBJDIR)/, DCDAnalyzer.o) -c $(addprefix $(DCDDIR)/, DCDAnalyzer.cpp)

$(addprefix $(OBJDIR)/, DCDReader.o): $(addprefix $(DCDDIR)/, DCDReader.cpp)
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix $(OBJDIR)/, DCDReader.o) -c $(addprefix $(DCDDIR)/, DCDReader.cpp)

$(addprefix $(OBJDIR)/, DCDWriter.o): $(addprefix $(DCDDIR)/, DCDWriter.cpp)
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix $(OBJDIR)/, DCDWriter.o) -c $(addprefix $(DCDDIR)/, DCDWriter.cpp)

$(addprefix $(OBJDIR)/, PDBParser.o): $(addprefix $(PDBDIR)/, PDBParser.cpp)
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix $(OBJDIR)/, PDBParser.o) -c $(addprefix $(PDBDIR)/, PDBParser.cpp)

$(addprefix $(OBJDIR)/, PDBReader.o): $(addprefix $(PDBDIR)/, PDBReader.cpp)
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix $(OBJDIR)/, PDBReader.o) -c $(addprefix $(PDBDIR)/, PDBReader.cpp)

$(addprefix $(OBJDIR)/, PDBWriter.o): $(addprefix $(PDBDIR)/, PDBWriter.cpp)
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix $(OBJDIR)/, PDBWriter.o) -c $(addprefix $(PDBDIR)/, PDBWriter.cpp)

$(addprefix $(OBJDIR)/, NinfoParser.o): $(addprefix $(NINFODIR)/, NinfoParser.cpp)
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix $(OBJDIR)/, NinfoParser.o) -c $(addprefix $(NINFODIR)/, NinfoParser.cpp)

$(addprefix $(OBJDIR)/, NinfoReader.o): $(addprefix $(NINFODIR)/, NinfoReader.cpp)
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix $(OBJDIR)/, NinfoReader.o) -c $(addprefix $(NINFODIR)/, NinfoReader.cpp)

$(addprefix $(OBJDIR)/, PSFParser.o): $(addprefix $(PSFDIR)/, PSFParser.cpp)
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix $(OBJDIR)/, PSFParser.o) -c $(addprefix $(PSFDIR)/, PSFParser.cpp)

$(addprefix $(OBJDIR)/, PSFReader.o): $(addprefix $(PSFDIR)/, PSFReader.cpp)
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix $(OBJDIR)/, PSFReader.o) -c $(addprefix $(PSFDIR)/, PSFReader.cpp)

#$(addprefix $(OBJDIR)/, FileOpenJudge.o): FileOpenJudge.cpp
#	$(CC) $(CFLAGS) $(INCLUDE) -c $(addprefix $(MISCDIR)/, FileOpenJudge.cpp)
$(addprefix $(OBJDIR)/, FileOpenJudge.o): $(addprefix $(MISCDIR)/, FileOpenJudge.cpp)
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix $(OBJDIR)/, FileOpenJudge.o) -c $(addprefix $(MISCDIR)/, FileOpenJudge.cpp)

$(addprefix $(OBJDIR)/, DensityFunction.o): $(addprefix $(MISCDIR)/, DensityFunction.cpp)
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix $(OBJDIR)/, DensityFunction.o) -c $(addprefix $(MISCDIR)/, DensityFunction.cpp)

# mains

$(addprefix $(OBJDIR)/, slideAnalyze.o): slideAnalyze.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix $(OBJDIR)/, slideAnalyze.o) -c slideAnalyze.cpp

$(addprefix $(OBJDIR)/, makePDBfromDCD.o): makePDBfromDCD.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix $(OBJDIR)/, makePDBfromDCD.o) -c makePDBfromDCD.cpp

$(addprefix $(OBJDIR)/, searchNativeContact.o): searchNativeContact.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix $(OBJDIR)/, searchNativeContact.o) -c searchNativeContact.cpp

$(addprefix $(OBJDIR)/, nativeContactedRange.o): nativeContactedRange.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix $(OBJDIR)/, nativeContactedRange.o) -c nativeContactedRange.cpp

$(addprefix $(OBJDIR)/, convertHistData.o): convertHistData.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix $(OBJDIR)/, convertHistData.o) -c convertHistData.cpp

$(addprefix $(OBJDIR)/, calcQscore.o): calcQscore.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix $(OBJDIR)/, calcQscore.o) -c calcQscore.cpp

$(addprefix $(OBJDIR)/, searchContactedHistoneChain.o): searchContactedHistoneChain.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix $(OBJDIR)/, searchContactedHistoneChain.o) -c searchContactedHistoneChain.cpp

$(addprefix $(OBJDIR)/, getProteinResiContacted2DNA.o): getProteinResiContacted2DNA.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix $(OBJDIR)/, getProteinResiContacted2DNA.o) -c getProteinResiContacted2DNA.cpp

$(addprefix $(OBJDIR)/, searchUnwrappedExtent.o): searchUnwrappedExtent.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix $(OBJDIR)/, searchUnwrappedExtent.o) -c searchUnwrappedExtent.cpp

$(addprefix $(OBJDIR)/, searchDNAContacted2Dyad.o): searchDNAContacted2Dyad.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix $(OBJDIR)/, searchDNAContacted2Dyad.o) -c searchDNAContacted2Dyad.cpp

$(addprefix $(OBJDIR)/, helicaseProcessingExtent.o): helicaseProcessingExtent.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix $(OBJDIR)/, helicaseProcessingExtent.o) -c helicaseProcessingExtent.cpp

$(addprefix $(OBJDIR)/, calcPDNSQscore.o): calcPDNSQscore.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix $(OBJDIR)/, calcPDNSQscore.o) -c calcPDNSQscore.cpp

$(addprefix $(OBJDIR)/, sortPDNSResidue.o): sortPDNSResidue.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix $(OBJDIR)/, sortPDNSResidue.o) -c sortPDNSResidue.cpp

$(addprefix $(OBJDIR)/, detectDNAProteinContacts.o): detectDNAProteinContacts.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix $(OBJDIR)/, detectDNAProteinContacts.o) -c detectDNAProteinContacts.cpp

$(addprefix $(OBJDIR)/, generateDNAProteinHeatMap.o): generateDNAProteinHeatMap.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix $(OBJDIR)/, generateDNAProteinHeatMap.o) -c generateDNAProteinHeatMap.cpp

$(addprefix $(OBJDIR)/, calcContactStateProbability.o): calcContactStateProbability.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix $(OBJDIR)/, calcContactStateProbability.o) -c calcContactStateProbability.cpp

$(addprefix $(OBJDIR)/, getPDNS_ELEContacts.o): getPDNS_ELEContacts.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix $(OBJDIR)/, getPDNS_ELEContacts.o) -c getPDNS_ELEContacts.cpp

$(addprefix $(OBJDIR)/, showStateTransition.o): showStateTransition.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix $(OBJDIR)/, showStateTransition.o) -c showStateTransition.cpp

$(addprefix $(OBJDIR)/, test.o): test.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix $(OBJDIR)/, test.o) -c test.cpp

