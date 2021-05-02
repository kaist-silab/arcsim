build/debug/util.o build/release/util.o: src/util.cpp src/util.hpp src/mesh.hpp src/transformation.hpp \
  src/spline.hpp src/vectors.hpp src/winport.hpp src/io.hpp
src/util.cpp src/util.hpp src/mesh.hpp src/transformation.hpp :
  src/spline.hpp src/vectors.hpp src/winport.hpp src/io.hpp :
