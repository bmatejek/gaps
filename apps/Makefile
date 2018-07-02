#
# Makefile for GAPS apps
#

opt:
	$(MAKE) target "TARGET=$@"

_debug:
	$(MAKE) target "TARGET=$@"

clean:
	$(MAKE) target "TARGET=$@"

release:
	mkdir -p ../release/apps
	cp Makefile ../release/apps
	$(MAKE) target "TARGET=$@"

target: 
	cd 3dview; $(MAKE) $(TARGET)
	cd downsampleview; $(MAKE) $(TARGET)
	cd ebrogold; $(MAKE) $(TARGET)
	cd ebroview; $(MAKE) $(TARGET)
	cd goldview; $(MAKE) $(TARGET)
	cd h5view; $(MAKE) $(TARGET)
	cd neuronview; $(MAKE) $(TARGET)
	cd nuclearview; $(MAKE) $(TARGET)
	cd overlay; $(MAKE) $(TARGET)
	cd segmentlists; $(MAKE) $(TARGET)
	cd skeletonview; $(MAKE) $(TARGET)
	cd skeletonvisualizer; $(MAKE) $(TARGET)
	cd skelmulticut; $(MAKE) $(TARGET)
	cd skelteaser; $(MAKE) $(TARGET)
#	cd skelview; $(MAKE) $(TARGET)
