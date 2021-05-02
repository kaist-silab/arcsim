build/debug/spline.o build/release/spline.o: src/spline.cpp src/spline.hpp src/vectors.hpp src/winport.hpp \
  src/util.hpp src/mesh.hpp src/transformation.hpp
src/spline.cpp src/spline.hpp src/vectors.hpp src/winport.hpp :
  src/util.hpp src/mesh.hpp src/transformation.hpp :
