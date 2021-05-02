build/debug/geometry.o build/release/geometry.o: src/geometry.cpp src/geometry.hpp src/mesh.hpp \
  src/transformation.hpp src/spline.hpp src/vectors.hpp src/winport.hpp \
  src/util.hpp
src/geometry.cpp src/geometry.hpp src/mesh.hpp :
  src/transformation.hpp src/spline.hpp src/vectors.hpp src/winport.hpp :
  src/util.hpp :
