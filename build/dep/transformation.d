build/debug/transformation.o build/release/transformation.o: src/transformation.cpp src/transformation.hpp \
  src/spline.hpp src/vectors.hpp src/winport.hpp
src/transformation.cpp src/transformation.hpp :
  src/spline.hpp src/vectors.hpp src/winport.hpp :
