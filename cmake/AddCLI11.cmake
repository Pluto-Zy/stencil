if (NOT DEFINED CLI11_PATH)
    set(CLI11_PATH "${CMAKE_SOURCE_DIR}/CLI11")
endif()

message(STATUS "Looking for CLI11 sources in ${CLI11_PATH}")

if (EXISTS "${CLI11_PATH}" AND IS_DIRECTORY "${CLI11_PATH}" AND EXISTS "${CLI11_PATH}/CMakeLists.txt")
    message(STATUS "Looking for CLI11 sources - Found")
    add_subdirectory(${CLI11_PATH})
else()
    message(STATUS "Looking for CLI11 sources - Not Found")
    message(STATUS "Did not find CLI11 sources. Fetching from GitHub...")

    include(FetchContent)
    FetchContent_Declare(
        CLI11
        # version 2.4.1
        URL https://github.com/CLIUtils/CLI11/archive/refs/tags/v2.4.1.zip
        URL_HASH MD5=96aadbda330b7a84ec91c42eec43aba4
    )
    FetchContent_MakeAvailable(CLI11)
endif()
