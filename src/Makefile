CC=g++
CFLAG= -std=c++17 -pedantic -Wall -g

DCDDIR= ../include/DCD
PDBDIR= ../include/PDB
NINFODIR= ../include/NINFO
PSFDIR= ../include/PSF
MISCDIR= ../include/misc
EXEDIR= ../bin

INCLUDE= -I. -I../include -I../include/misc

DCDANALYSIS= DCDAnalyzer.o DCDParser.o
DCDREAD= DCDReader.o DCDParser.o
PDBREAD= PDBReader.o PDBParser.o
PDBWRITE= PDBWriter.o PDBParser.o
NINFOREAD= NinfoReader.o NinfoParser.o
PSFREAD= PSFReader.o PSFParser.o
FILEJUDGE= FileOpenJudge.o
MAKEDENSITY= DensityFunction.o
SLIDEMAIN= slideAnalyze.o 
PDBFROMDCD= makePDBfromDCD.o
NATIVECONTACT= searchNativeContact.o
NCRANGE= nativeContactedRange.o
CONVERTHIST= convertHistData.o
CALCQSCORE= calcQScore.o
CONTACTHISTONE= searchContactedHistoneChain.o
GETPDNS= getProteinResiContacted2DNA.o
SEARCHUNWRAP= searchUnwrappedExtent.o
SEARCHDNA2DYAD= searchDNAContacted2Dyad.o
HELICASEPROCESSX= helicaseProcessingExtent.o
CALCPDNSQSCORE= calcPDNSQscore.o
SORTPDNSRESIDUE= sortPDNSResidue.o
DETECTDNAPROTEINCONTACTS= detectDNAProteinContacts.o
GENERATEDNAPROHEATMAP= generateDNAProteinHeatMap.o
CALCCONTACTSTATEPROB= calcContactStateProbability.o
GETPDNSELECON= getPDNS_ELEContacts.o

TEST= test/read_psf.o
ALL= ../include/DCD/DCDAnalyzer.o ../include/DCD/DCDParser.o ../include/DCD/DCDReader.o ../include/PDB/PDBParser.o ../include/PDB/PDBReader.o ../include/PDB/PDBWriter.o ../include/misc/FileOpenJudge.o ../include/misc/DensityFunction.o slideAnalyze.o makePDBfromDCD.o searchNativeContact.o convertHistData.o nativeContactedRange.o ../include/NINFO/NinfoParser.o ../include/NINFO/NinfoReader.o calcQScore.o searchContactedHistoneChain.o getProteinResiContacted2DNA.o searchUnwrappedExtent.o searchDNAContacted2Dyad.o helicaseProcessingExtent.o calcPDNSQscore.o sortPDNSResidue.o detectDNAProteinContacts.o generateDNAProteinHeatMap.o calcContactStateProbability.o getPDNS_ELEContacts.o

SLIDEEXE= slide
PDBFROMDCDEXE= dcd2pdb
SEARCHNATCONEXE= searchnative
SEARCHHISTONECONTACTEXE= searchcontacthistone
NATCONRANGEEXE= ncrange
HISTGRAMEXE= histgram
QSCOREEXE= qscore
GETPDNSEXE= getpdns
SEARCHUNWRAPEXE= unwrapextent
SLIDE2DYADEXE= slide2dyad
HELICASEPROCESSXEXE= helicaseprocess
CALCPDNSQSCOREEXE= pdnsqscore
SORTPDNSRESIDUEEXE= sortpdns
DETECTDNAPROTEINCONTACTSEXE= detectdnaprocontact
GENERATEDNAPROHEATMAPEXE= pdnsmapping
CALCCONTACTSTATEPROBEXE= contactstate2prob
GETPDNSELECONEXE= getcontactaddress

TESTEXE= test

SLIDEANALYZE= $(SLIDEMAIN) $(addprefix $(MISCDIR)/, $(FILEJUDGE)) $(addprefix $(DCDDIR)/, $(DCDANALYSIS)) $(addprefix $(PDBDIR)/, $(PDBREAD))

DCD2PDB= $(PDBFROMDCD) $(addprefix $(MISCDIR)/, $(FILEJUDGE)) $(addprefix $(DCDDIR)/, $(DCDREAD)) $(addprefix $(PDBDIR)/, $(PDBWRITE))

SEARCHNATCON= $(NATIVECONTACT) $(addprefix $(MISCDIR)/, $(FILEJUDGE)) $(addprefix $(DCDDIR)/, $(DCDANALYSIS)) $(addprefix $(PDBDIR)/, $(PDBREAD))

NATCONRANGE= $(NCRANGE) $(addprefix $(MISCDIR)/, $(FILEJUDGE)) $(addprefix $(DCDDIR)/, $(DCDANALYSIS)) $(addprefix $(PDBDIR)/, $(PDBREAD))

HISTGRAM= $(addprefix $(MISCDIR)/, $(MAKEDENSITY)) $(CONVERTHIST)

TESTSET= $(TEST) $(addprefix $(PSFDIR)/, $(PSFREAD))

QSCORE= $(CALCQSCORE) $(addprefix $(MISCDIR)/, $(FILEJUDGE)) $(addprefix $(NINFODIR)/, $(NINFOREAD)) $(addprefix $(DCDDIR)/, $(DCDANALYSIS))

SEARCHHISCON= $(CONTACTHISTONE) $(addprefix $(MISCDIR)/, $(FILEJUDGE)) $(addprefix $(DCDDIR)/, $(DCDANALYSIS)) $(addprefix $(PSFDIR)/, $(PSFREAD))

GETPDNSRESI= $(GETPDNS) $(addprefix $(MISCDIR)/, $(FILEJUDGE)) $(addprefix $(NINFODIR)/, $(NINFOREAD)) $(addprefix $(PSFDIR)/, $(PSFREAD)) $(addprefix $(PDBDIR)/, $(PDBREAD))

SEARCHUNWRAPEXTENT= $(SEARCHUNWRAP) $(addprefix $(MISCDIR)/, $(FILEJUDGE)) $(addprefix $(DCDDIR)/, $(DCDANALYSIS)) $(addprefix $(PSFDIR)/, $(PSFREAD)) $(addprefix $(NINFODIR)/, $(NINFOREAD))

SLIDE2DYAD= $(SEARCHDNA2DYAD) $(addprefix $(DCDDIR)/, $(DCDANALYSIS)) $(addprefix $(PSFDIR)/, $(PSFREAD))

HELICASEPROCESS= $(HELICASEPROCESSX) $(addprefix $(DCDDIR)/, $(DCDANALYSIS)) $(addprefix $(PSFDIR)/, $(PSFREAD))

PDNSQSCORE= $(CALCPDNSQSCORE) $(addprefix $(NINFODIR)/, $(NINFOREAD)) $(addprefix $(DCDDIR)/, $(DCDANALYSIS)) $(addprefix $(PSFDIR)/, $(PSFREAD))

SORTPDNS= $(SORTPDNSRESIDUE) $(addprefix $(DCDDIR)/, $(DCDANALYSIS)) $(addprefix $(PSFDIR)/, $(PSFREAD)) $(addprefix $(NINFODIR)/, $(NINFOREAD))

#DETECTDNAPROCONTACT= $(DETECTDNAPROTEINCONTACTS) $(addprefix $(NINFODIR)/, $(NINFOREAD)) $(addprefix $(DCDDIR)/, $(DCDANALYSIS)) $(addprefix $(PSFDIR)/, $(PSFREAD))
DETECTDNAPROCONTACT= $(DETECTDNAPROTEINCONTACTS) $(addprefix $(DCDDIR)/, $(DCDANALYSIS)) $(addprefix $(PSFDIR)/, $(PSFREAD))

PDNSMAPPING= $(GENERATEDNAPROHEATMAP) $(addprefix $(DCDDIR)/, $(DCDANALYSIS)) $(addprefix $(PSFDIR)/, $(PSFREAD))

CONTACTSTATE2PROB= $(CALCCONTACTSTATEPROB) $(addprefix $(PSFDIR)/, $(PSFREAD))

GETCONTACTADDRESS= $(GETPDNSELECON) $(addprefix $(PSFDIR)/, $(PSFREAD)) $(addprefix $(PDBDIR)/, $(PDBREAD))


slide: $(SLIDEANALYZE)
	$(CC) $(INCLUDE) -o $(addprefix $(EXEDIR)/, $(SLIDEEXE)) $(SLIDEANALYZE)

dcd2pdb: $(DCD2PDB)
	$(CC) $(INCLUDE) -o $(addprefix $(EXEDIR)/, $(PDBFROMDCDEXE)) $(DCD2PDB)

searchnative: $(SEARCHNATCON)
	$(CC) $(INCLUDE) -o $(addprefix $(EXEDIR)/, $(SEARCHNATCONEXE)) $(SEARCHNATCON)

ncrange: $(NATCONRANGE)
	$(CC) $(INCLUDE) -o $(addprefix $(EXEDIR)/, $(NATCONRANGEEXE)) $(NATCONRANGE)

histgram: $(HISTGRAM)
	$(CC) $(INCLUDE) -o $(addprefix $(EXEDIR)/, $(HISTGRAMEXE)) $(HISTGRAM)

test: $(TESTSET)
	$(CC) $(INCLUDE) -o $(addprefix $(EXEDIR)/, $(TESTEXE)) $(TESTSET)

qscore: $(QSCORE)
	$(CC) $(INCLUDE) -o $(addprefix $(EXEDIR)/, $(QSCOREEXE)) $(QSCORE)

searchcontacthistone: $(SEARCHHISCON)
	$(CC) $(INCLUDE) -o $(addprefix $(EXEDIR)/, $(SEARCHHISTONECONTACTEXE)) $(SEARCHHISCON)

getpdns: $(GETPDNSRESI)
	$(CC) $(INCLUDE) -o $(addprefix $(EXEDIR)/, $(GETPDNSEXE)) $(GETPDNSRESI)

unwrapextent: $(SEARCHUNWRAPEXTENT)
	$(CC) $(INCLUDE) -o $(addprefix $(EXEDIR)/, $(SEARCHUNWRAPEXE)) $(SEARCHUNWRAPEXTENT)

slide2dyad: $(SLIDE2DYAD)
	$(CC) $(INCLUDE) -o $(addprefix $(EXEDIR)/, $(SLIDE2DYADEXE)) $(SLIDE2DYAD)

helicaseprocess: $(HELICASEPROCESS)
	$(CC) $(INCLUDE) -o $(addprefix $(EXEDIR)/, $(HELICASEPROCESSXEXE)) $(HELICASEPROCESS)

pdnsqscore: $(PDNSQSCORE)
	$(CC) $(INCLUDE) -o $(addprefix $(EXEDIR)/, $(CALCPDNSQSCOREEXE)) $(PDNSQSCORE)

sortpdns: $(SORTPDNS)
	$(CC) $(INCLUDE) -o $(addprefix $(EXEDIR)/, $(SORTPDNSRESIDUEEXE)) $(SORTPDNS)


detectdnaprocontact: $(DETECTDNAPROCONTACT)
	$(CC) $(INCLUDE) -o $(addprefix $(EXEDIR)/, $(DETECTDNAPROTEINCONTACTSEXE)) $(DETECTDNAPROCONTACT)

pdnsmapping: $(PDNSMAPPING)
	$(CC) $(INCLUDE) -o $(addprefix $(EXEDIR)/, $(GENERATEDNAPROHEATMAPEXE)) $(PDNSMAPPING)

contactstate2prob: $(CONTACTSTATE2PROB)
	$(CC) $(INCLUDE) -o $(addprefix $(EXEDIR)/, $(CALCCONTACTSTATEPROBEXE)) $(CONTACTSTATE2PROB)

getcontactaddress: $(GETCONTACTADDRESS)
	$(CC) $(INCLUDE) -o $(addprefix $(EXEDIR)/, $(GETPDNSELECONEXE)) $(GETCONTACTADDRESS)

clean:
	rm $(ALL)
