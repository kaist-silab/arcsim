build/debug/misc.o build/release/misc.o: src/misc.cpp src/io.hpp src/mesh.hpp src/transformation.hpp \
  src/spline.hpp src/vectors.hpp src/winport.hpp src/util.hpp
src/misc.cpp src/io.hpp src/mesh.hpp src/transformation.hpp :
  src/spline.hpp src/vectors.hpp src/winport.hpp src/util.hpp :
