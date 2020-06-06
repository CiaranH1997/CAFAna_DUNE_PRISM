set(scripts_to_build
  #  PRISMPred_SystProf_ChiSqSurface
  #  PRISMPred_SystProf_MySurface
  MakePRISMPredInterps
  OffAxisEventRates
  OffAxisNDCAFCombiner
  PRISMPrediction
  PRISMDChi2Scan
  DumpEventRateTable
  #  ChiSqSurfacePRISMPred
  PRISMRegOptimizer
  PRISMPred_ChiSqSurf
)

foreach(SCRPT ${scripts_to_build})
  install(FILES ${SCRPT}.C DESTINATION scripts)

  add_executable(${SCRPT} ${SCRPT}.C)

  target_link_libraries(${SCRPT} ${CAFANALIBS})

    target_link_libraries(${SCRPT} ${CAFANA_LINK_DIRS})
    target_link_libraries(${SCRPT} ${CAFANA_DEPEND_LIBS})

  if(DEFINED USE_GPERFTOOLS AND USE_GPERFTOOLS)
    add_dependencies(${SCRPT} gperftools)
    target_link_libraries(${SCRPT} ${GPERF_LINKER_FLAGS})
  endif()

  if(DEFINED USE_OPENMP AND USE_OPENMP)
    target_compile_options(${SCRPT} BEFORE PUBLIC -DUSE_PREDINTERP_OMP -fopenmp)
  endif()

  install(TARGETS ${SCRPT} DESTINATION bin)
endforeach()
