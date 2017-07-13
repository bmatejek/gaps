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
	cd ebrogold; $(MAKE) $(TARGET)
	cd featureview; $(MAKE) $(TARGET)
	cd h5view; $(MAKE) $(TARGET)
	cd neuronview; $(MAKE) $(TARGET)
	cd skel_dist; $(MAKE) $(TARGET)
	cd skeleton; $(MAKE) $(TARGET)
	cd skeletonview; $(MAKE) $(TARGET)
