cmake_minimum_required(VERSION 3.15)

# conditional adding of subdirectories
function(custom_add_subdir dirname)
    if (IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${dirname})
        message("adding ${dirname}")
        add_subdirectory(${dirname})
    else()
        message("skip ${dirname}")
    endif()
endfunction(custom_add_subdir dirname)

#custom_add_subdir(PartToGridConverter)
#custom_add_subdir(FuncAtPoint)
#custom_add_subdir(WaterProfile)
#custom_add_subdir(SPH2D_Drawer)
#custom_add_subdir(SPH2D_PicGen)

if (WIN32)
    set(SPH2DParamsGeneratorPath ${CMAKE_SOURCE_DIR}/tools/SPH2DParamsGenerator/SPH2DParamsGenerator/bin/Release/)
    if (EXISTS ${SPH2DParamsGeneratorPath})
        #install(DIRECTORY "${SPH2DParamsGeneratorPath}" DESTINATION "RRSPH/SPH2DParamsGenerator")
    endif()
endif()

install(DIRECTORY "${CMAKE_SOURCE_DIR}/tools/Scripts" DESTINATION "RRSPH")