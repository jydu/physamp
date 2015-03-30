#! /bin/sh
arch=`uname -m`
version=0.2.0-1

strip physamp/bppalnoptim
tar cvzf physamp-${arch}-bin-static-${version}.tar.gz physamp/bppalnoptim

