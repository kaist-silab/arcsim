build/debug/dde.o build/release/dde.o: src/dde.cpp src/dde.hpp src/sparse.hpp src/util.hpp src/mesh.hpp \
  src/transformation.hpp src/spline.hpp src/vectors.hpp src/winport.hpp \
  src/cloth.hpp
src/dde.cpp src/dde.hpp src/sparse.hpp src/util.hpp src/mesh.hpp :
  src/transformation.hpp src/spline.hpp src/vectors.hpp src/winport.hpp :
  src/cloth.hpp :
