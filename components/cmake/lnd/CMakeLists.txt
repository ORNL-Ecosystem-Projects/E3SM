build_model(${COMP_CLASS} ${COMP_NAME})

option(BUILD_ELM_OFFLINE_DRIVER "Build the standalone ELM offline driver executable" OFF)

if (BUILD_ELM_OFFLINE_DRIVER)
  add_executable(elm_offline_driver
    ${SRCROOT}/components/elm/src/cpl/elm_offline_driver.F90)

  set_target_properties(elm_offline_driver PROPERTIES LINKER_LANGUAGE CXX)
  target_link_libraries(elm_offline_driver PRIVATE lnd csm_share)
  if (CMAKE_CXX_IMPLICIT_LINK_LIBRARIES)
    target_link_libraries(elm_offline_driver PRIVATE ${CMAKE_CXX_IMPLICIT_LINK_LIBRARIES})
  endif()
  if (USE_PETSC)
    target_link_libraries(elm_offline_driver PRIVATE "${PETSC_LIBRARIES}")
    target_include_directories(elm_offline_driver PRIVATE "${PETSC_INCLUDES}")
  endif()
  if (NOT DEFINED ENV{SKIP_BLAS})
    target_link_libraries(elm_offline_driver PRIVATE BLAS::BLAS LAPACK::LAPACK)
  endif()
  target_include_directories(elm_offline_driver PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}
    ${EXEROOT}/cmake-bld/cmake/lnd
    ${EXEROOT}/cmake-bld/cmake/csm_share)
endif()
