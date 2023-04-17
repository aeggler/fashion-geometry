# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/buildNew/_deps/glad-src"
  "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/buildNew/_deps/glad-build"
  "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/buildNew/_deps/glad-subbuild/glad-populate-prefix"
  "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/buildNew/_deps/glad-subbuild/glad-populate-prefix/tmp"
  "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/buildNew/_deps/glad-subbuild/glad-populate-prefix/src/glad-populate-stamp"
  "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/buildNew/_deps/glad-subbuild/glad-populate-prefix/src"
  "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/buildNew/_deps/glad-subbuild/glad-populate-prefix/src/glad-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/buildNew/_deps/glad-subbuild/glad-populate-prefix/src/glad-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/Users/annaeggler/Desktop/Masterarbeit/fashion-descriptors/buildNew/_deps/glad-subbuild/glad-populate-prefix/src/glad-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
