<!--
SPDX-FileCopyrightText: 2026 Julien Y. Dutheil <dutheil@evolbio.mpg.de>

SPDX-License-Identifier: GPL-3.0-or-later
-->

optimize the alignment using various methods:

```bash
bppalnoptim param=AlnOptimFasttree.bpp
bppalnoptim param=AlnOptimAuto.bpp LINKAGE=complete
bppalnoptim param=AlnOptimAuto.bpp LINKAGE=single
bppalnoptim param=AlnOptimAuto.bpp LINKAGE=average
bppalnoptim param=AlnOptimAuto.bpp LINKAGE=median
bppalnoptim param=AlnOptimAuto.bpp LINKAGE=centroid
bppalnoptim param=AlnOptimAuto.bpp LINKAGE=ward
```

Get a selection of sequences:
```bash
bppalnoptim param=AlnOptimFasttreeStop.bpp
```
