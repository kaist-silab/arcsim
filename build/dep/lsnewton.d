build/debug/lsnewton.o build/release/lsnewton.o: src/lsnewton.cpp src/optimization.hpp src/sparse.hpp \
  src/util.hpp src/mesh.hpp src/transformation.hpp src/spline.hpp \
  src/vectors.hpp src/winport.hpp src/taucs.hpp
src/lsnewton.cpp src/optimization.hpp src/sparse.hpp :
  src/util.hpp src/mesh.hpp src/transformation.hpp src/spline.hpp :
  src/vectors.hpp src/winport.hpp src/taucs.hpp :
