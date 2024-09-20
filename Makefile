.PHONY: all 
all: libdocopt.a energyAnalysis partialGridEnergy

libdocopt.a:
	$(MAKE) -C docopt.cpp

.PHONY: energyAnalysis
energyAnalysis: libdocopt.a
	$(MAKE) -C energyAnalysis

.PHONY: partialGridEnergy
partialGridEnergy: libdocopt.a
	$(MAKE) -C partialGridEnergy

.PHONY: clean
clean:
	$(MAKE) -C docopt.cpp clean
	$(MAKE) -C energyAnalysis clean
	$(MAKE) -C partialGridEnergy clean



