CPPFLAGS=`root-config --cflags` -I/usr/local/root/include
LDFLAGS=`root-config --libs` -lMinuit
#anaZgTree: anaZgTree.cc
skimmer: skimmer.cc
	g++ $(CPPFLAGS) $(LDFLAGS) -o $@ $^

