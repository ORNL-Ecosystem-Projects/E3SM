add_subdirectory(base)
add_subdirectory(betr)
add_subdirectory(fates)
add_subdirectory(stub)
add_subdirectory(alquimia)

set(EMI_LIBRARIES emi_em_base; emi_em_betr; emi_em_fates; emi_em_alquimia; emi_em_stub; ${EMI_LIBRARIES} PARENT_SCOPE)

