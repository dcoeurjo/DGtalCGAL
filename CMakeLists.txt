#David Coeurjolly <david.coeurjolly@liris.cnrs.fr>
# -----------------------------------------------------------------------------
# Project definition
# -----------------------------------------------------------------------------
PROJECT(DGtalCGAL)

cmake_minimum_required (VERSION 2.6) 


set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)
if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()


find_package(CGAL COMPONENTS Core )
include( ${CGAL_USE_FILE} )


FIND_PACKAGE(DGtal 0.5 REQUIRED)
INCLUDE_DIRECTORIES(${DGTAL_INCLUDE_DIRS})
LINK_DIRECTORIES(${DGTAL_LIBRARY_DIRS})

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR})



FILE(GLOB_RECURSE DGTAL_MYHEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*h)
add_custom_target(_headers SOURCES ${DGTAL_MYHEADERS})

SET(DGTAL_TOOLS_SRC
  implicitToPCL
  #curvatureFromImplicit
)


FOREACH(FILE ${DGTAL_TOOLS_SRC})
  add_executable(${FILE} ${FILE})
  target_link_libraries (${FILE}  ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${DGTAL_LIBRARIES})
  install(TARGETS ${FILE} 
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
ENDFOREACH(FILE)

