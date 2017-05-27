# CMake script for physamp package
# Author: Julien Dutheil
# Created: 24/10/2014

ADD_EXECUTABLE(bppalnoptim AlignmentOptimizer.cpp)
ADD_EXECUTABLE(bppphysamp PhyloSampler.cpp)
SET (physamp-targets
	bppalnoptim
	bppphysamp
	)

FOREACH (target ${physamp-targets})
  # Link (static or shared)
  IF (BUILD_STATIC)
    TARGET_LINK_LIBRARIES (${target} ${BPP_LIBS_STATIC})
    SET_TARGET_PROPERTIES (${target} LINK_SEARCH_END_STATIC TRUE)
  ELSE (BUILD_STATIC)
    TARGET_LINK_LIBRARIES (${target} ${BPP_LIBS_SHARED})
  ENDIF (BUILD_STATIC)
ENDFOREACH (target)

INSTALL (TARGETS ${physamp-targets} DESTINATION bin)
