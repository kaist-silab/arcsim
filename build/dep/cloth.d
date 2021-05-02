build/debug/cloth.o build/release/cloth.o: src/cloth.cpp src/cloth.hpp src/dde.hpp src/sparse.hpp \
  src/util.hpp src/mesh.hpp src/transformation.hpp src/spline.hpp \
  src/vectors.hpp src/winport.hpp
src/cloth.cpp src/cloth.hpp src/dde.hpp src/sparse.hpp :
  src/util.hpp src/mesh.hpp src/transformation.hpp src/spline.hpp :
  src/vectors.hpp src/winport.hpp :
