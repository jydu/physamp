# SPDX-FileCopyrightText: Julien Y. Dutheil
#
# SPDX-License-Identifier: CECILL-2.1

add_executable(bppalnoptim AlignmentOptimizer.cpp)
add_executable(bppphysamp PhyloSampler.cpp)
set(physamp-targets bppalnoptim bppphysamp)

foreach(target ${physamp-targets})
    # Link (static or shared)
    if(BUILD_STATIC)
        target_link_libraries(${target} ${BPP_LIBS_STATIC})
        set_target_properties(${target} PROPERTIES LINK_SEARCH_END_STATIC TRUE)
	target_link_options(${target} PRIVATE -static-libgcc -static-libstdc++)
    else(BUILD_STATIC)
        target_link_libraries(${target} ${BPP_LIBS_SHARED})
        set_target_properties(
            ${target}
            PROPERTIES
                INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}"
        )
    endif(BUILD_STATIC)
endforeach(target)

install(TARGETS ${physamp-targets} DESTINATION ${CMAKE_INSTALL_BINDIR})
