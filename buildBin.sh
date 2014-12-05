#! /bin/sh
arch=x86_64 #i686
version=0.1.0-1

strip physamp/bppalnoptim
tar cvzf physamp-${arch}-bin-static-${version}.tar.gz physamp/bppalnoptim

