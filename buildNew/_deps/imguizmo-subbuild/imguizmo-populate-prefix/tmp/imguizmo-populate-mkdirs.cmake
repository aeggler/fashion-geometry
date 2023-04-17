# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/buildNew/_deps/imguizmo-src"
  "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/buildNew/_deps/imguizmo-build"
  "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/buildNew/_deps/imguizmo-subbuild/imguizmo-populate-prefix"
  "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/buildNew/_deps/imguizmo-subbuild/imguizmo-populate-prefix/tmp"
  "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/buildNew/_deps/imguizmo-subbuild/imguizmo-populate-prefix/src/imguizmo-populate-stamp"
  "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/buildNew/_deps/imguizmo-subbuild/imguizmo-populate-prefix/src"
  "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/buildNew/_deps/imguizmo-subbuild/imguizmo-populate-prefix/src/imguizmo-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/buildNew/_deps/imguizmo-subbuild/imguizmo-populate-prefix/src/imguizmo-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/buildNew/_deps/imguizmo-subbuild/imguizmo-populate-prefix/src/imguizmo-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
