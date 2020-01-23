all:
	rootcint -f Dict_datamonitor.cpp -c `root-config --cflags` -I${ONLINE_MAIN}/include -I/home/eic/analysis/mpgd4eic/mapping/include -p DataMonitor.h LinkDef.h
	g++ -o bin/datamonitor gui.cc DataMonitor.cxx Dict_datamonitor.cpp -I${ONLINE_MAIN}/include -I/home/eic/analysis/mpgd4eic/mapping/include -L${ONLINE_MAIN}/lib -L/home/eic/analysis/mpgd4eic/lib/lib -lmapping -lpmonitor -lEvent -lNoRootEvent -lmessage `root-config --cflags --glibs` -fpermissive
	rm Dict_datamonitor.cpp
	cp Dict_datamonitor_rdict.pcm bin
