###INCLUDE_BEGIN### /cvmfs/larsoft.opensciencegrid.org/products/cetmodules/v3_21_01/config/package-config.cmake.preamble.in (2 passes)
########################################################################
# sbnanaConfig.cmake
#
#   Config for CMake project sbnana 09.78.06.
#
# Generated by cetmodules 3.21.01 at Tue Aug 20 17:18:25 CDT 2024.
#
# Compiled from:
#   /cvmfs/larsoft.opensciencegrid.org/products/cetmodules/v3_21_01/config/package-config.cmake.preamble.in
#   /cvmfs/larsoft.opensciencegrid.org/products/cetmodules/v3_21_01/config/package-config.cmake.init.in
#   genConfig/sbnana-generated-config.cmake.vars.in
#   genConfig/sbnana-generated-config.cmake.deps.in
#   genConfig/sbnana-generated-config.cmake.targets.in
#   genConfig/sbnana-generated-config.cmake.target_vars.in
#   /cvmfs/larsoft.opensciencegrid.org/products/cetmodules/v3_21_01/config/package-config.cmake.bottom.in
########################################################################
###INCLUDE_END###   /cvmfs/larsoft.opensciencegrid.org/products/cetmodules/v3_21_01/config/package-config.cmake.preamble.in

###INCLUDE_BEGIN### /cvmfs/larsoft.opensciencegrid.org/products/cetmodules/v3_21_01/config/package-config.cmake.init.in (2 passes)
###############################################################################
####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was sbnanaConfig.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../../" ABSOLUTE)

macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

macro(check_required_components _NAME)
  foreach(comp ${${_NAME}_FIND_COMPONENTS})
    if(NOT ${_NAME}_${comp}_FOUND)
      if(${_NAME}_FIND_REQUIRED_${comp})
        set(${_NAME}_FOUND FALSE)
      endif()
    endif()
  endforeach()
endmacro()

####################################################################################

# Make find_dependency() macro available.
include(CMakeFindDependencyMacro)
###INCLUDE_END###   /cvmfs/larsoft.opensciencegrid.org/products/cetmodules/v3_21_01/config/package-config.cmake.init.in

###INCLUDE_BEGIN### genConfig/sbnana-generated-config.cmake.vars.in
####################################
# Package variable definitions.
####################################
# sbnana_NAMESPACE
set(sbnana_NAMESPACE "sbnana")
# sbnana_INCLUDE_DIR
if (EXISTS "${PACKAGE_PREFIX_DIR}/include")
  file(GLOB _sbnana_TMP_DIR_ENTRIES "${PACKAGE_PREFIX_DIR}/include/*")
  if (_sbnana_TMP_DIR_ENTRIES)
    set_and_check(sbnana_INCLUDE_DIR "${PACKAGE_PREFIX_DIR}/include")
  endif()
  unset(_sbnana_TMP_DIR_ENTRIES)
endif()
# sbnana_LIBRARY_DIR
if (EXISTS "${PACKAGE_PREFIX_DIR}/slf7.x86_64.e20.prof/lib")
  file(GLOB _sbnana_TMP_DIR_ENTRIES "${PACKAGE_PREFIX_DIR}/slf7.x86_64.e20.prof/lib/*")
  if (_sbnana_TMP_DIR_ENTRIES)
    set_and_check(sbnana_LIBRARY_DIR "${PACKAGE_PREFIX_DIR}/slf7.x86_64.e20.prof/lib")
  endif()
  unset(_sbnana_TMP_DIR_ENTRIES)
endif()
# sbnana_BIN_DIR
if (EXISTS "${PACKAGE_PREFIX_DIR}/slf7.x86_64.e20.prof/bin")
  file(GLOB _sbnana_TMP_DIR_ENTRIES "${PACKAGE_PREFIX_DIR}/slf7.x86_64.e20.prof/bin/*")
  if (_sbnana_TMP_DIR_ENTRIES)
    set_and_check(sbnana_BIN_DIR "${PACKAGE_PREFIX_DIR}/slf7.x86_64.e20.prof/bin")
  endif()
  unset(_sbnana_TMP_DIR_ENTRIES)
endif()
# sbnana_INSTALLED_SOURCE_DIR
if (EXISTS "${PACKAGE_PREFIX_DIR}/source")
  file(GLOB _sbnana_TMP_DIR_ENTRIES "${PACKAGE_PREFIX_DIR}/source/*")
  if (_sbnana_TMP_DIR_ENTRIES)
    set_and_check(sbnana_INSTALLED_SOURCE_DIR "${PACKAGE_PREFIX_DIR}/source")
  endif()
  unset(_sbnana_TMP_DIR_ENTRIES)
endif()

####################################
# Package include directories.
####################################
if (IS_DIRECTORY "${PACKAGE_PREFIX_DIR}/include")
  # CMake convention:
  set(sbnana_INCLUDE_DIRS "${PACKAGE_PREFIX_DIR}/include")
endif()

####################################
# Package library directories.
####################################
if (IS_DIRECTORY "${PACKAGE_PREFIX_DIR}/slf7.x86_64.e20.prof/lib")
  # CMake convention:
  set(sbnana_LIBRARY_DIRS "${PACKAGE_PREFIX_DIR}/slf7.x86_64.e20.prof/lib")
endif()
###INCLUDE_END###   genConfig/sbnana-generated-config.cmake.vars.in

###INCLUDE_BEGIN### genConfig/sbnana-generated-config.cmake.deps.in
####################################
# Transitive dependencies.
####################################
set(_sbnana_PACKAGE_PREFIX_DIR "${PACKAGE_PREFIX_DIR}")
cet_find_library(IFDH NAMES ifdh PATHS ENV IFDHC_LIB NO_DEFAULT_PATH)
set(PACKAGE_PREFIX_DIR "${_sbnana_PACKAGE_PREFIX_DIR}")
unset(_sbnana_PACKAGE_PREFIX_DIR)
###INCLUDE_END###   genConfig/sbnana-generated-config.cmake.deps.in

###INCLUDE_BEGIN### genConfig/sbnana-generated-config.cmake.targets.in
####################################
# Exported targets, and package components.
####################################

##################
# Automatically-generated runtime targets: sbnana
##################
include("${CMAKE_CURRENT_LIST_DIR}/sbnanaTargets.cmake")
foreach (component IN LISTS sbnana_FIND_COMPONENTS)
  include("${CMAKE_CURRENT_LIST_DIR}/sbnanaTargets_${component}.cmake")
endforeach()

##################
# Manually-generated non-runtime targets.
##################
set(_targetsDefined)
set(_targetsNotDefined)
set(_expectedTargets)
foreach (_expectedTarget IN ITEMS
    sbnana::cafana.py
    sbnana::rootlogon.C
    sbnana::cafe
    sbnana::load_cafana_libs.C)
  list(APPEND _expectedTargets ${_expectedTarget})
  if (NOT TARGET ${_expectedTarget})
    list(APPEND _targetsNotDefined ${_expectedTarget})
  endif()
  if (TARGET ${_expectedTarget})
    list(APPEND _targetsDefined ${_expectedTarget})
  endif()
endforeach()
if ("${_targetsDefined}" STREQUAL "${_expectedTargets}")  # Nothing to do.
elseif ("${_targetsDefined}" STREQUAL "") # Need to define targets.
  add_executable(sbnana::cafana.py IMPORTED)
  set_target_properties(sbnana::cafana.py
    PROPERTIES IMPORTED_LOCATION "${PACKAGE_PREFIX_DIR}/slf7.x86_64.e20.prof/bin/cafana.py")
  add_executable(sbnana::rootlogon.C IMPORTED)
  set_target_properties(sbnana::rootlogon.C
    PROPERTIES IMPORTED_LOCATION "${PACKAGE_PREFIX_DIR}/slf7.x86_64.e20.prof/bin/rootlogon.C")
  add_executable(sbnana::cafe IMPORTED)
  set_target_properties(sbnana::cafe
    PROPERTIES IMPORTED_LOCATION "${PACKAGE_PREFIX_DIR}/slf7.x86_64.e20.prof/bin/cafe")
  add_executable(sbnana::load_cafana_libs.C IMPORTED)
  set_target_properties(sbnana::load_cafana_libs.C
    PROPERTIES IMPORTED_LOCATION "${PACKAGE_PREFIX_DIR}/slf7.x86_64.e20.prof/bin/load_cafana_libs.C")
else()
  message(FATAL_ERROR "Some (but not all) targets in this export set were already defined.
Targets Defined: ${_targetsDefined}
Targets not yet defined: ${_targetsNotDefined}
")
endif()
unset(_targetsDefined)
unset(_targetsNotDefined)
unset(_expectedTargets)
###INCLUDE_END###   genConfig/sbnana-generated-config.cmake.targets.in

###INCLUDE_BEGIN### genConfig/sbnana-generated-config.cmake.target_vars.in
####################################
# Old cetbuildtools-style target variables.
####################################
if (${CETMODULES_CURRENT_PROJECT_NAME}_OLD_STYLE_CONFIG_VARS OR # Per-dependent setting.
 cetbuildtools_UPS_VERSION OR # Backward-compatibility.
 cetbuildtools IN_LIST ${CETMODULES_CURRENT_PROJECT_NAME}_UPS_BUILD_ONLY_DEPENDENCIES)
  set(SBNANA_CAFANACORE sbnana::CAFAnaCore)
  set(SBNANA_CAFANAUNFOLD sbnana::CAFAnaUnfold)
  set(SBNANA_CAFANAVARS sbnana::CAFAnaVars)
  set(SBNANA_CAFANACUTS sbnana::CAFAnaCuts)
  set(SBNANA_CAFANASYSTS sbnana::CAFAnaSysts)
  set(SBNANA_CAFANAEXTRAP sbnana::CAFAnaExtrap)
  set(SBNANA_CAFANAPREDICTION sbnana::CAFAnaPrediction)
  set(SBNANA_CAFANAEXPERIMENT sbnana::CAFAnaExperiment)
  set(SBNANA_CAFANAXSEC sbnana::CAFAnaXSec)
  set(SBNANA_CAFANAANALYSIS sbnana::CAFAnaAnalysis)
  set(SBNANA_SBNANAVARS sbnana::SBNAnaVars)
  set(SBNANA_SBNANACUTS sbnana::SBNAnaCuts)
endif()
###INCLUDE_END###   genConfig/sbnana-generated-config.cmake.target_vars.in

###INCLUDE_BEGIN### /cvmfs/larsoft.opensciencegrid.org/products/cetmodules/v3_21_01/config/package-config.cmake.bottom.in
# Check package components.
check_required_components(sbnana)
###INCLUDE_END###   /cvmfs/larsoft.opensciencegrid.org/products/cetmodules/v3_21_01/config/package-config.cmake.bottom.in