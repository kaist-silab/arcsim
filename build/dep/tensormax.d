build/debug/tensormax.o build/release/tensormax.o: src/tensormax.cpp src/tensormax.hpp src/vectors.hpp \
  src/winport.hpp src/util.hpp src/mesh.hpp src/transformation.hpp \
  src/spline.hpp
src/tensormax.cpp src/tensormax.hpp src/vectors.hpp :
  src/winport.hpp src/util.hpp src/mesh.hpp src/transformation.hpp :
  src/spline.hpp :
