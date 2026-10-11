remove_definitions(-DRDKIT_MOLDRAW2D_BUILD)

rdkit_nanobind_extension(rdMolDraw2D
                         rdMolDraw2D.cpp
                         DEST Chem/Draw
                         LINK_LIBRARIES MolDraw2D )

add_pytest(pyMolDraw2D ${CMAKE_CURRENT_SOURCE_DIR}/../Wrap/testMolDraw2D.py)
add_pytest(pyMechanismDrawing ${CMAKE_CURRENT_SOURCE_DIR}/../Wrap/testMechanismDrawing.py)
add_pytest(pyMechanismAutomatic ${CMAKE_CURRENT_SOURCE_DIR}/../Wrap/testMechanismAutomatic.py)
add_pytest(pyMechanismGallery ${CMAKE_CURRENT_SOURCE_DIR}/../Wrap/testMechanismGallery.py)
