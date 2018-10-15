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
	cd 2dview; $(MAKE) $(TARGET)
	cd 3dview; $(MAKE) $(TARGET)
	cd downsampleview; $(MAKE) $(TARGET)
	cd edgeview; $(MAKE) $(TARGET)
	cd errorview; $(MAKE) $(TARGET)
	cd featureview; $(MAKE) $(TARGET)
	cd nodeview; $(MAKE) $(TARGET)
	cd overlay; $(MAKE) $(TARGET)
	cd segmentlists; $(MAKE) $(TARGET)
	cd skeletonbenchmark; $(MAKE) $(TARGET)
	cd skeletonedit; $(MAKE) $(TARGET)
	cd skeletonview; $(MAKE) $(TARGET)
