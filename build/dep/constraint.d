build/debug/constraint.o build/release/constraint.o: src/constraint.cpp src/constraint.hpp src/mesh.hpp \
  src/transformation.hpp src/spline.hpp src/vectors.hpp src/winport.hpp \
  src/sparse.hpp src/util.hpp src/magic.hpp
src/constraint.cpp src/constraint.hpp src/mesh.hpp :
  src/transformation.hpp src/spline.hpp src/vectors.hpp src/winport.hpp :
  src/sparse.hpp src/util.hpp src/magic.hpp :
