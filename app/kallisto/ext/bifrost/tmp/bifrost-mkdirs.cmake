# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION ${CMAKE_VERSION}) # this file comes with cmake

# If CMAKE_DISABLE_SOURCE_CHANGES is set to true and the source directory is an
# existing directory in our source tree, calling file(MAKE_DIRECTORY) on it
# would cause a fatal error, even though it would be a no-op.
if(NOT EXISTS "/home/adriantkachenko/Programming/Research/FASTR/app/kallisto/ext/bifrost")
  file(MAKE_DIRECTORY "/home/adriantkachenko/Programming/Research/FASTR/app/kallisto/ext/bifrost")
endif()
file(MAKE_DIRECTORY
  "/home/adriantkachenko/Programming/Research/FASTR/app/kallisto/ext/bifrost/src/bifrost-build"
  "/home/adriantkachenko/Programming/Research/FASTR/app/kallisto/ext/bifrost"
  "/home/adriantkachenko/Programming/Research/FASTR/app/kallisto/ext/bifrost/tmp"
  "/home/adriantkachenko/Programming/Research/FASTR/app/kallisto/ext/bifrost/src/bifrost-stamp"
  "/home/adriantkachenko/Programming/Research/FASTR/app/kallisto/ext/bifrost/src"
  "/home/adriantkachenko/Programming/Research/FASTR/app/kallisto/ext/bifrost/src/bifrost-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/home/adriantkachenko/Programming/Research/FASTR/app/kallisto/ext/bifrost/src/bifrost-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/home/adriantkachenko/Programming/Research/FASTR/app/kallisto/ext/bifrost/src/bifrost-stamp${cfgdir}") # cfgdir has leading slash
endif()
