build/debug/bvh.o build/release/bvh.o: src/bvh.cpp src/bvh.hpp src/mesh.hpp src/transformation.hpp \
  src/spline.hpp src/vectors.hpp src/winport.hpp src/util.hpp \
  src/collision.hpp src/cloth.hpp src/dde.hpp src/sparse.hpp \
  src/constraint.hpp
src/bvh.cpp src/bvh.hpp src/mesh.hpp src/transformation.hpp :
  src/spline.hpp src/vectors.hpp src/winport.hpp src/util.hpp :
  src/collision.hpp src/cloth.hpp src/dde.hpp src/sparse.hpp :
  src/constraint.hpp :
