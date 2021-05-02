build/debug/vectors.o build/release/vectors.o: src/vectors.cpp src/vectors.hpp src/winport.hpp \
  src/blockvectors.hpp src/util.hpp src/mesh.hpp src/transformation.hpp \
  src/spline.hpp
src/vectors.cpp src/vectors.hpp src/winport.hpp :
  src/blockvectors.hpp src/util.hpp src/mesh.hpp src/transformation.hpp :
  src/spline.hpp :
