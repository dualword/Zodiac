# CMakeLists.txt for the main sourcecode of molsKetch

# Including qt4 and OpenBabel
set(QT_USE_QTASSISTANT TRUE)
set(QT_USE_QTSVG TRUE)
include(${QT_USE_FILE})
include_directories(${OPENBABEL2_INCLUDE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

# Set the variable with the .ui files
set(molsketch_UIS settings.ui)
qt4_wrap_ui(molsketch_UIS_H ${molsketch_UIS})

# Set the variable with all the sourcecode
set(molsketch_SRCS main.cpp molecule.cpp atom.cpp bond.cpp element.cpp mollibitem.cpp mainwindow.cpp molview.cpp molscene.cpp settings.cpp fileio.cpp commands.cpp periodictablewidget.cpp)

# Including the resources
qt4_add_resources(molsketch_SRCS molsketch.qrc)

# Add a program icon for windows
if(WIN32)
  set(molsketch_SRCS ${molsketch_SRCS} icon.rc)
endif(WIN32)

# MOC-ing headers
set(molsketch_MOC_HDRS mainwindow.h molscene.h molview.h settings.h periodictablewidget.h)
qt4_wrap_cpp(molsketch_SRCS ${molsketch_MOC_HDRS})

# Take care of the .moc files
qt4_automoc(${molsketch_SRCS})

# Create an executable from this source
add_executable(molsketch ${molsketch_SRCS} ${molsketch_UIS_H})

# Link the code against Qt and OpenBabel
target_link_libraries(molsketch ${QT_LIBRARIES} ${OPENBABEL2_LIBRARIES})


# Install the executable
install(TARGETS molsketch DESTINATION bin)

# Install the library
install(DIRECTORY ${PROJECT_SOURCE_DIR}/library DESTINATION share/molsketch PATTERN ".svn" EXCLUDE)

# Install the documentation
install(DIRECTORY ${PROJECT_SOURCE_DIR}/doc DESTINATION share/doc/molsketch PATTERN ".svn" EXCLUDE)

# Install menu entries on Linux
if(UNIX)
  install(FILES ${PROJECT_SOURCE_DIR}/src/molsketch.desktop DESTINATION share/applications)
  install(FILES ${PROJECT_SOURCE_DIR}/src/images/molsketch.xpm DESTINATION share/pixmaps)
  install(FILES ${PROJECT_SOURCE_DIR}/src/images/molsketch.png DESTINATION share/icons/hicolor/128x128/apps)
endif(UNIX)

# Check subdirs
add_subdirectory(i18n)
if(KDE4_FOUND)
  add_subdirectory(part)
endif(KDE4_FOUND)