# add to sieve subfolder
add_subdirectory ("dgs")

# p2 heuristic (only for test purpose)
add_executable (test_p2_heuristic p2_heuristic.cpp)
install(TARGETS test_p2_heuristic RUNTIME DESTINATION ${BINSUFFIX}/sieve)


# test sampler (only for test purpose)
set (SAMPLER_SOURCES sampler_basic.cpp)
add_executable(test_sampler
  sampler_main.cpp
  dgs/dgs_gauss_dp.c
  dgs/dgs_gauss_mp.c
  dgs/dgs_bern.c
  ${SAMPLER_SOURCES}
  ${Arith_Sources}
  )
target_link_libraries (test_sampler -lgmp -lmpfr -lqd)
install(TARGETS test_sampler RUNTIME DESTINATION ${BINSUFFIX}/sieve)

# Gauss sieve (MP09)
set (SIEVE_GAUSS_SOURCES sieve_Gauss.cpp)

add_executable(test_sieve
  sieve_main.cpp
  ${SIEVE_GAUSS_SOURCES}
  dgs/dgs_gauss_dp.c
  dgs/dgs_gauss_mp.c
  dgs/dgs_bern.c
  ${SAMPLER_SOURCES}
  ${Arith_Sources}
)
target_link_libraries (test_sieve -lgmp -lmpfr -lqd)
install(TARGETS test_sieve RUNTIME DESTINATION ${BINSUFFIX}/sieve)

# 3-sieve
