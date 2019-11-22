CC=g++
CFLAG= -std=c++14 -pedantic -Wall -g
OBJSLIDE= slide_analyze.o DCDAnalyzer.o DCDParser.o
OBJLENGTH= length.o DCDAnalyzer.o DCDParser.o 
OBJNATCONTACT= analyze_native_contact.o DCDAnalyzer.o DCDParser.o FileOpenJudge.o

slideanalyze: $(OBJSLIDE)
	$(CC) -o slideanalyze $(OBJSLIDE)

lenanalyze: $(OBJLENGTH)
	$(CC) -o length $(OBJLENGTH)

natcontact: $(OBJNATCONTACT)
	$(CC) -o natcontact $(OBJNATCONTACT)

clean:
	rm DCDAnalyzer.o DCDParser.o FileOpenJudge.o
