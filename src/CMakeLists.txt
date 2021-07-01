set (SOURCES
     Hydroinfo_h5.cpp
     Arsenal.cpp
     ThermalPhoton.cpp
     Table2D.cpp
     PhotonEmission.cpp
     ParameterReader.cpp
     tensor_trans.cpp
     gauss_quadrature.cpp
     Hydroinfo_MUSIC.cpp
     HadronGasRhoSpectralFunction.cpp
    )

add_library(${libname} SHARED ${SOURCES})
if (link_with_lib)
    add_executable(${exename} main.cpp)
else (link_with_lib)
    add_executable(${exename} main.cpp ${SOURCES})
endif (link_with_lib)
target_link_libraries(${libname} ${HDF5_LIBRARIES})

target_link_libraries(${exename} ${libname} ${HDF5_LIBRARIES})

install(TARGETS ${exename} ${libname}
        RUNTIME DESTINATION ${CMAKE_HOME_DIRECTORY}
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib/static)
