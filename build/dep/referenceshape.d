build/debug/referenceshape.o build/release/referenceshape.o: src/referenceshape.cpp src/referenceshape.hpp \
  src/mesh.hpp src/transformation.hpp src/spline.hpp src/vectors.hpp \
  src/winport.hpp
src/referenceshape.cpp src/referenceshape.hpp :
  src/mesh.hpp src/transformation.hpp src/spline.hpp src/vectors.hpp :
  src/winport.hpp :
