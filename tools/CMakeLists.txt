cmake_minimum_required(VERSION 3.15)

add_subdirectory(PartToGridConverter)
add_subdirectory(SPH2D_FIO2)
add_subdirectory(FuncAtPoint)
add_subdirectory(WaterProfile)
add_subdirectory(SPH2D_Drawer)
add_subdirectory(SPH2D_PicGen)


set(SPH2DParamsGeneratorPath ${CMAKE_SOURCE_DIR}/tools/SPH2DParamsGenerator)
if (EXISTS ${SPH2DParamsGeneratorPath})
    install(DIRECTORY "${SPH2DParamsGeneratorPath}/SPH2DParamsGenerator/bin/Release/" DESTINATION "SPH2D/SPH2DParamsGenerator")
endif()