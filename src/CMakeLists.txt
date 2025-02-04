set (SOURCES
     Hydroinfo_h5.cpp
     Arsenal.cpp
     ThermalPhoton.cpp
     ThermalDilepton.cpp
     Table2D.cpp
     PhotonEmission.cpp
     ParameterReader.cpp
     tensor_trans.cpp
     gauss_quadrature.cpp
     Hydroinfo_MUSIC.cpp
     QGP2to2Total.cpp
     QGPAMYCollinear.cpp
     HadronGasRhoSpectralFunction.cpp
     HadronGasPipiBremsstrahlung.cpp
     HadronGasPiRhoOmega.cpp
     DileptonQGPLO.cpp
    )

add_library(${libname} SHARED ${SOURCES})
target_link_libraries(${libname} ${HDF5_LIBRARIES})
if (link_with_lib)
    add_executable(${exename} main.cpp)
    target_link_libraries(${exename} ${libname} ${HDF5_LIBRARIES})
else (link_with_lib)
    add_executable(${exename} main.cpp ${SOURCES})
    target_link_libraries(${exename} ${HDF5_LIBRARIES})
endif (link_with_lib)

install(TARGETS ${exename} ${libname}
        RUNTIME DESTINATION ${CMAKE_HOME_DIRECTORY}
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib/static)
