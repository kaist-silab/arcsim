#!/bin/bash
# Installer script for ARCSim 0.3.1 on Ubuntu and other Debian-based distributions
RED='\033[0;31m'
GREEN='\033[0;32m'
ORANGE='\033[0;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${GREEN}Would you like to install the required packages for Ubuntu / Debian-based distributions?${NC}"
select yn in "Yes" "No"
do
	case $yn in
		Yes) echo "Installing necessary packages..."; apt-get update; 
			apt-get install gcc make g++ libboost-all-dev freeglut3-dev \
				gfortran liblapacke-dev libpng-dev libpng++-dev scons \
        			libatlas-base-dev ctags libopenblas-dev; break;;
		No) echo "Skipping..."; break;;
	esac
done

cd dependencies/
echo "Running make on dependencies folder. Grab a cup of coffee, this will take a while."
make
if [ $? -ne 0 ]
then
	echo -e "${RED}Make error on dependencies occurred. Exiting...${NC}"
	exit 1
fi

cd ../
echo "Modifying line 118 in src/sparse.hpp..."
sed -i '118s|.*|    file << "}]" << std::endl;|' src/sparse.hpp
echo "Commenting line 805 and 811 in dependencies/include/taucs.h..."
sed -i '805s|.*|/*extern int isnan(double);*/|' dependencies/include/taucs.h
sed -i '811s|.*|/*extern int isinf(double);*/|' dependencies/include/taucs.h

echo "Running make on main folder. This will also take a while."
make -f "Makefile.linux"
if [ $? -ne 0 ]
then
	echo -e "${RED}Make error on main folder occurred. Exiting...${NC}"
	exit 1
fi

echo -e "${GREEN}stallation finished. Would you like to test it?${NC}"
select yn in "Yes" "No"
do
	case $yn in
		Yes ) echo "Testing flag simulation..."; 
			./bin/arcsim simulate conf/flag.json; break;;
		No ) echo "Skipping..."; break;;
	esac
done
echo "Exiting script..."

