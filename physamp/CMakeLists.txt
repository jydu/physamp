# CMake script for physamp package
# Author: Julien Dutheil
# Created: 24/10/2014

ADD_EXECUTABLE(bppalnoptim AlignmentOptimizer.cpp)
TARGET_LINK_LIBRARIES(bppalnoptim ${LIBS})
SET_TARGET_PROPERTIES(bppalnoptim PROPERTIES LINK_SEARCH_END_STATIC ${BUILD_STATIC})

ADD_EXECUTABLE(bppphysamp PhyloSampler.cpp)
TARGET_LINK_LIBRARIES(bppphysamp ${LIBS})
SET_TARGET_PROPERTIES(bppphysamp PROPERTIES LINK_SEARCH_END_STATIC ${BUILD_STATIC})

# Install progs
INSTALL(TARGETS
  bppalnoptim
  bppphysamp
  DESTINATION bin)
