build/debug/mot_parser.o build/release/mot_parser.o: src/mot_parser.cpp src/mot_parser.hpp src/obstacle.hpp \
  src/mesh.hpp src/transformation.hpp src/spline.hpp src/vectors.hpp \
  src/winport.hpp src/util.hpp src/io.hpp
src/mot_parser.cpp src/mot_parser.hpp src/obstacle.hpp :
  src/mesh.hpp src/transformation.hpp src/spline.hpp src/vectors.hpp :
  src/winport.hpp src/util.hpp src/io.hpp :
