build/debug/morph.o build/release/morph.o: src/morph.cpp src/morph.hpp src/mesh.hpp src/transformation.hpp \
  src/spline.hpp src/vectors.hpp src/winport.hpp src/geometry.hpp \
  src/util.hpp
src/morph.cpp src/morph.hpp src/mesh.hpp src/transformation.hpp :
  src/spline.hpp src/vectors.hpp src/winport.hpp src/geometry.hpp :
  src/util.hpp :
