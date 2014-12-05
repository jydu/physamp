#! /bin/sh
arch=x86_64 #i686
version=0.1.0-1

strip bppSuite/bppalnoptim
tar cvzf physamp-${arch}-bin-static-${version}.tar.gz PhySamp/bppalnoptim

