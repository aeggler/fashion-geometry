# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/Users/annaeggler/Desktop/Documents/Studium/Masterarbeit/fashion-descriptors/build-new/_deps/eigen-src"
  "/Users/annaeggler/Desktop/Documents/Studium/Masterarbeit/fashion-descriptors/build-new/_deps/eigen-build"
  "/Users/annaeggler/Desktop/Documents/Studium/Masterarbeit/fashion-descriptors/build-new/_deps/eigen-subbuild/eigen-populate-prefix"
  "/Users/annaeggler/Desktop/Documents/Studium/Masterarbeit/fashion-descriptors/build-new/_deps/eigen-subbuild/eigen-populate-prefix/tmp"
  "/Users/annaeggler/Desktop/Documents/Studium/Masterarbeit/fashion-descriptors/build-new/_deps/eigen-subbuild/eigen-populate-prefix/src/eigen-populate-stamp"
  "/Users/annaeggler/Desktop/Documents/Studium/Masterarbeit/fashion-descriptors/build-new/_deps/eigen-subbuild/eigen-populate-prefix/src"
  "/Users/annaeggler/Desktop/Documents/Studium/Masterarbeit/fashion-descriptors/build-new/_deps/eigen-subbuild/eigen-populate-prefix/src/eigen-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/Users/annaeggler/Desktop/Documents/Studium/Masterarbeit/fashion-descriptors/build-new/_deps/eigen-subbuild/eigen-populate-prefix/src/eigen-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/Users/annaeggler/Desktop/Documents/Studium/Masterarbeit/fashion-descriptors/build-new/_deps/eigen-subbuild/eigen-populate-prefix/src/eigen-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
