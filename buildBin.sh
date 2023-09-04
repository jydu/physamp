#! /bin/sh
arch=`uname -m`
version=1.1.0-1

strip physamp/bppalnoptim
strip physamp/bppphysamp
tar cvzf physamp-${arch}-bin-static-${version}.tar.gz physamp/bppalnoptim physamp/bppphysamp

