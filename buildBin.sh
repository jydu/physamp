#! /bin/sh

# SPDX-FileCopyrightText: 2026 Julien Y. Dutheil <dutheil@evolbio.mpg.de>
#
# SPDX-License-Identifier: GPL-3.0-or-later

arch=`uname -m`
version=1.2.0-1

strip physamp/bppalnoptim
strip physamp/bppphysamp
tar cvzf physamp-${arch}-bin-static-${version}.tar.gz physamp/bppalnoptim physamp/bppphysamp

