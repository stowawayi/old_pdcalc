cmake_minimum_required(VERSION 3.10)
project(pdcalc Fortran)

# Set Fortran standard and flags
set(CMAKE_Fortran_STANDARD 77)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)

# Default to empty string
set(compiler_flags "")

# Compiler-specific flags
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(compiler_flags "-Wall -Wunused-parameter -fno-range-check -std=legacy -ffpe-trap=zero,denormal,invalid,overflow -C")

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    set(compiler_flags "-warn all -stand legacy -fpe0 -check bounds,uninit")

    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
        set(compiler_flags "${compiler_flags} -traceback")
    endif()

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Absoft")
    set(compiler_flags "-C -s -N113")

    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
        set(compiler_flags "${compiler_flags} -B")
    endif()  
endif()

# Apply the compiler flags
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${compiler_flags}")

# Source files for the static library
set(LIB_SOURCES
    acon.f
    dypres.f
    errmsg.f
    etcalc.f
    intgf.f
    lncalc.f
    ogh.f
    overp.f
    pdcalc.f
    ranf.f
    wrcalc.f
    wrclcy.f
    wrcrtr.f
    wrpers.f
)

# Create the static library
add_library(pd STATIC ${LIB_SOURCES})

# Main program source and driver sources
set(MAIN_SOURCE main.for)

# Add the executable and link it against the library
add_executable(main ${MAIN_SOURCE})
target_link_libraries(main PRIVATE pd)

# Optional install
install(TARGETS main DESTINATION bin)
install(TARGETS pd ARCHIVE DESTINATION lib)
